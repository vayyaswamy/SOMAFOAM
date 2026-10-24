// Mode diagnostics at the end of the step (createModeDiagnostics.H), 1D.
//
// The electron net frequency nu_net = d(electron source)/dNe is the sum over
// reactions of (electrons produced - consumed) times the rate per electron:
// ionization counts positive, attachment and recombination negative. It does
// not depend on Ne itself, so cells at the density floor give the right
// alpha_eff. lnM is the largest sum of alpha_eff dy over a contiguous interval
// (any start, any end), i.e. the largest gain an avalanche can reach crossing
// that interval in the present field.
if (modeDiagnostics1D)
{
    const label eI = modeDiagElectronI;

    const scalarField& Ne = mspm().N(eI).internalField();
    const scalarField mue(mspm().mu(eI).internalField());
    const scalarField nu(chemistry.dRRDi(eI)().internalField());
    const vectorField& Ec = E.internalField();

    const scalarField& magSfG = mesh.magSf().boundaryField()[groundPatchI];
    const scalar area = gSum(magSfG);
    const vector axis = gSum(mesh.Sf().boundaryField()[groundPatchI])/area;

    const scalar yD =
        gAverage(mesh.Cf().boundaryField()[modeDiagDielectricPatchI] & axis);
    const scalar yG = gAverage(mesh.Cf().boundaryField()[groundPatchI] & axis);
    const scalar gap = yG - yD;

    // Cells in order along the axis
    const scalarField yc((mesh.C().internalField() & axis) - yD);
    labelList order;
    sortedOrder(yc, order);

    const scalarField dy(mesh.V().field()/area);

    // Positive ions
    scalarField Nplus(Ne.size(), 0.0);
    forAll(mspm().species(), i)
    {
        if (mspm().z(i) > 0)
        {
            Nplus += mspm().N(i).internalField();
        }
    }

    scalar NeMax = 0, yNeMax = 0, NpMax = 0, EMax = 0, yEMax = 0;
    scalar sigmaP = 0, sigmaE = 0;
    scalar best = 0, run = 0, yBestEnd = 0;

    forAll(order, k)
    {
        const label c = order[k];
        const scalar magE = mag(Ec[c]);

        if (Ne[c] > NeMax)
        {
            NeMax = Ne[c];
            yNeMax = yc[c];
        }
        if (Nplus[c] > NpMax)
        {
            NpMax = Nplus[c];
        }
        if (magE > EMax)
        {
            EMax = magE;
            yEMax = yc[c];
        }

        sigmaP += Nplus[c]*dy[c];
        sigmaE += Ne[c]*dy[c];

        // Largest interval sum of alpha_eff dy (Kadane)
        const scalar a = nu[c]/max(mue[c]*magE, VSMALL)*dy[c];
        run = max(run + a, a);
        if (run > best)
        {
            best = run;
            yBestEnd = yc[c];
        }
    }

    const scalar e = plasmaConstants::eCharge.value();
    sigmaP *= e;
    sigmaE *= e;

    // Wall fields, positive towards the gas, and the mean gap field
    const scalarField& magSfD =
        mesh.magSf().boundaryField()[modeDiagDielectricPatchI];
    const scalar Vg =
        gSum(Phi.boundaryField()[groundPatchI]*magSfG)/area;
    const scalar Vd =
        gSum(Phi.boundaryField()[modeDiagDielectricPatchI]*magSfD)
       /gSum(magSfD);

    const scalar EGround =
        -gSum((E.boundaryField()[groundPatchI] & axis)*magSfG)/area;
    const scalar EDiel =
        gSum((E.boundaryField()[modeDiagDielectricPatchI] & axis)*magSfD)
       /gSum(magSfD);

    if (modeDiagFile.valid())
    {
        modeDiagFile()
            << runTime.value() << tab << runTime.deltaTValue() << tab
            << NeMax << tab << yNeMax/gap << tab
            << NpMax << tab << EMax << tab << yEMax/gap << tab
            << mag(Vg - Vd)/gap << tab << EGround << tab << EDiel << tab
            << best << tab << yBestEnd/gap << tab
            << sigmaP << tab << sigmaE << endl;
    }
}
