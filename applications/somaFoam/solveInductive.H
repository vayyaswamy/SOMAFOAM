// Solve the azimuthal vector potential with the current plasma conductivity and update
// the cycle-averaged inductive power density, scaled to the prescribed absorbed power.
if (isICP && (icpPending || icpVerify || runTime.timeIndex() % icpInterval == 0))
{
    const fvMesh& emMesh = emMeshPtr();
    volScalarField& AR = ARPtr();
    volScalarField& AI = AIPtr();
    volVector2Field& A = APtr();
    const scalar w = 2*mathematicalConstant::pi*icpFreq;
    const scalar mu0 = 4e-7*mathematicalConstant::pi;
    const scalar qe = 1.602176634e-19, meM = 9.1093837e-31;

    // plasma cell -> em cell map, rebuilt if the plasma mesh changed
    {
        const vectorField& Cp = mesh.C().internalField();
        const vectorField& Cem = emMesh.C().internalField();
        if (emCellOfPlasma.size() != Cp.size())
        {
            emCellOfPlasma.setSize(Cp.size(), -1);
            label nSearch = 0;
            forAll(Cp, i)
            {
                if (i < Cem.size() && mag(Cem[i] - Cp[i]) < 1e-9) { emCellOfPlasma[i] = i; }
                else { emCellOfPlasma[i] = emMesh.findCell(Cp[i]); nSearch++; }
                if (emCellOfPlasma[i] < 0)
                {
                    FatalError << "ICP: plasma cell " << i << " at " << Cp[i]
                        << " is outside the em mesh" << abort(FatalError);
                }
            }
            Info<< "ICP: mapped " << Cp.size() << " plasma cells, "
                << nSearch << " by search" << endl;
        }
    }

    const scalarField& Ne = mspm().N(icpEIndex).internalField();
    const scalarField& muE = mspm().mu(icpEIndex).internalField();

    volScalarField sigR(IOobject("sigR", runTime.timeName(), emMesh), emMesh,
        dimensionedScalar("zero", dimless/dimArea, 0.0));
    volScalarField sigI(IOobject("sigI", runTime.timeName(), emMesh), emMesh,
        dimensionedScalar("zero", dimless/dimArea, 0.0));

    // The electron mobility is set by the first plasma solve; until then the
    // collision frequency is undefined and the field update is postponed
    const scalar muMin = returnReduce(muE.size() ? min(mag(muE)) : GREAT, minOp<scalar>());
    const bool ready = icpVerify || (muMin > SMALL);
    icpPending = !ready;

    if (icpVerify)
    {
        // analytic plasma n0 J0(2.405 r/R) sin(pi z/L) with uniform nuM, for comparison
        // with the independent solvers of the preliminary results
        const vectorField& Cem = emMesh.C().internalField();
        forAll(Cem, k)
        {
            const scalar r = Foam::sqrt(sqr(Cem[k].x()) + sqr(Cem[k].y())), z = Cem[k].z();
            if (r < icpVR && z < icpVL)
            {
                const scalar x = 2.405*r/icpVR;
                scalar j0 = 0, term = 1;
                for (label m = 0; m < 20; m++) { j0 += term; term *= -sqr(x/2)/sqr(scalar(m + 1)); }
                const scalar n = max(icpVn0*j0*Foam::sin(mathematicalConstant::pi*z/icpVL), 0.0);
                const scalar s0 = n*qe*qe/meM/(sqr(icpVnu) + sqr(w));
                sigR.internalField()[k] = w*mu0*s0*icpVnu;
                sigI.internalField()[k] = -w*mu0*s0*w;
            }
        }
    }
    else if (ready)
    {
        forAll(emCellOfPlasma, i)
        {
            const scalar nuM = qe/(meM*mag(muE[i]));
            const scalar s0 = max(Ne[i], 0.0)*qe*qe/meM/(sqr(nuM) + sqr(w));
            // scaled by w mu0 so that the equations below are dimensionally uniform
            sigR.internalField()[emCellOfPlasma[i]] = w*mu0*s0*nuM;
            sigI.internalField()[emCellOfPlasma[i]] = -w*mu0*s0*w;
        }
    }
    else
    {
        Info<< "ICP: electron mobility not yet set, field update postponed" << endl;
    }

    // decide whether the field must be re-solved. The measure is the change of the
    // complex conductivity since the last solve weighted by the field intensity of that
    // solve, which is the first-order relative change of the plasma current and power.
    const label stepNow = runTime.timeIndex();
    scalar dSigma = GREAT;
    if (ready && !icpVerify && icpLastSolve >= 0 && icpSigRLast.size() == sigR.size())
    {
        const scalarField& Vem = emMesh.V();
        scalar num = 0, den = 0;
        forAll(sigR, k)
        {
            const scalar a2V = (sqr(AR.internalField()[k]) + sqr(AI.internalField()[k]))*Vem[k];
            num += Foam::sqrt(sqr(sigR.internalField()[k] - icpSigRLast[k])
                + sqr(sigI.internalField()[k] - icpSigILast[k]))*a2V;
            den += Foam::sqrt(sqr(icpSigRLast[k]) + sqr(icpSigILast[k]))*a2V;
        }
        num = returnReduce(num, sumOp<scalar>());
        den = returnReduce(den, sumOp<scalar>());
        dSigma = (den > VSMALL) ? num/den : GREAT;
    }
    const bool solveNow = ready &&
    (
        icpVerify || icpTol <= 0 || icpLastSolve < 0
     || stepNow - icpLastSolve >= icpMaxInterval || dSigma > icpTol
    );

    if (ready && !solveNow)
    {
        // field kept; only the coil current follows the power ramp
        const scalar Pset = (icpRamp > 0)
            ? icpPower*min(1.0, runTime.value()/icpRamp) : icpPower;
        if (icpPsetLast > 0 && Pset != icpPsetLast)
        {
            scalarField& P = PindPtr().internalField();
            P *= Pset/icpPsetLast;
            icpPindPrev = P;
            PindPtr().correctBoundaryConditions();
            icpCoilI = (icpP1 > VSMALL) ? Foam::sqrt(Pset/icpP1) : 0;
        }
        icpPsetLast = Pset;
    }

    if (solveNow)
    {
    if (!icpVerify)
    {
        Info<< "ICP: field update at step " << stepNow << " after "
            << (icpLastSolve < 0 ? 0 : stepNow - icpLastSolve) << " steps, conductivity change "
            << (dSigma < GREAT ? dSigma : 0) << ", n_e " << gMin(Ne) << " to " << gMax(Ne)
            << " m^-3, nu_m " << qe/(meM*returnReduce(max(mag(muE)), maxOp<scalar>()))
            << " to " << qe/(meM*muMin) << " 1/s" << endl;
    }
    icpSigRLast = sigR.internalField();
    icpSigILast = sigI.internalField();
    icpLastSolve = stepNow;
    icpNSolves++;
    volScalarField invR2(IOobject("invR2", runTime.timeName(), emMesh), emMesh,
        dimensionedScalar("zero", dimless/dimArea, 0.0));
    {
        const vectorField& Cem = emMesh.C().internalField();
        forAll(Cem, i)
        {
            invR2.internalField()[i] = 1.0/max(sqr(Cem[i].x()) + sqr(Cem[i].y()), 1e-12);
        }
    }
    volScalarField src(IOobject("src", runTime.timeName(), emMesh), emMesh,
        dimensionedScalar("zero", AR.dimensions()/dimArea, 0.0));
    src.internalField() = mu0*JcoilEM;

    // (1/r) d/dr(r dA/dr) + d2A/dz2 - A/r^2 - i w mu0 sigma A = -mu0 J, block-coupled in (AR, AI)
    fvScalarMatrix eqR(fvm::laplacian(AR) - fvm::Sp(invR2, AR) + fvm::Sp(sigI, AR) + src);
    fvScalarMatrix eqI(fvm::laplacian(AI) - fvm::Sp(invR2, AI) + fvm::Sp(sigI, AI));
    fvScalarMatrix cRI(fvm::Sp(sigR, AI));
    fvScalarMatrix cIR(fvm::Sp(-sigR, AR));

    fvBlockMatrix<vector2> AEqn(A);
    AEqn.insertEquation(0, eqR);
    AEqn.insertEquation(1, eqI);
    AEqn.insertEquationCoupling(0, 1, cRI);
    AEqn.insertEquationCoupling(1, 0, cIR);
    BlockSolverPerformance<vector2> perf = AEqn.solve();
    AEqn.retrieveSolution(0, AR.internalField());
    AEqn.retrieveSolution(1, AI.internalField());
    AR.correctBoundaryConditions();
    AI.correctBoundaryConditions();

    // power per (1 A)^2 and scaling of the coil current to the prescribed power
    scalarField& P = PindPtr().internalField();
    const scalarField& V = mesh.V();
    scalar P1 = 0;
    forAll(emCellOfPlasma, i)
    {
        const label k = emCellOfPlasma[i];
        const scalar sR = sigR.internalField()[k]/(w*mu0);
        P[i] = 0.5*sR*sqr(w)*(sqr(AR.internalField()[k]) + sqr(AI.internalField()[k]));
        P1 += P[i]*V[i];
    }
    P1 = returnReduce(P1, sumOp<scalar>())*360.0/icpWedge;
    const scalar Pset = (icpRamp > 0)
        ? icpPower*min(1.0, runTime.value()/icpRamp) : icpPower;
    icpCoilI = (P1 > VSMALL) ? Foam::sqrt(Pset/P1) : 0;
    P *= sqr(icpCoilI);
    if (icpRelax < 1 && icpPindPrev.size() == P.size())
    {
        P = icpRelax*P + (1 - icpRelax)*icpPindPrev;
    }
    icpPindPrev = P;
    icpP1 = P1;
    icpPsetLast = Pset;
    PindPtr().correctBoundaryConditions();
    Info<< "ICP: solve " << icpNSolves << ", block solve iterations " << perf.nIterations()
        << ", set power " << Pset << " W, power per A^2 " << P1 << " W, coil current " << icpCoilI
        << " A, peak Pind " << gMax(P) << " W/m3" << endl;

    if (icpVerify)
    {
        // power per (1 A)^2 and peak |E| per 1 A in the plasma region from the em mesh
        scalar Pv = 0, Ev = 0;
        const vectorField& Cem = emMesh.C().internalField();
        const scalarField& Vem = emMesh.V();
        forAll(Cem, k)
        {
            const scalar a2 = sqr(AR.internalField()[k]) + sqr(AI.internalField()[k]);
            Pv += 0.5*sigR.internalField()[k]/(w*mu0)*sqr(w)*a2*Vem[k];
            if (Cem[k].z() < icpVL) { Ev = max(Ev, w*Foam::sqrt(a2)); }
        }
        Pv = returnReduce(Pv, sumOp<scalar>())*360.0/icpWedge;
        Ev = returnReduce(Ev, maxOp<scalar>());
        Info<< "ICP verification: power per (1 A)^2 " << Pv << " W, peak |E| per 1 A "
            << Ev << " V/m, block iterations " << perf.nIterations() << endl;
        AR.write(); AI.write();
        return 0;
    }
    }
}
