// Total current density at the end of the step, 1D.
//
// The Poisson equation of a step uses the charge of the previous step, and
// eps*dE/dt between steps amplifies small step-to-step changes of the field
// at steps of 1e-11 s, so it is noisy and lags the conduction current.
//
// In 1D the total current density j = Jc + eps dE/dt is uniform across gas
// and dielectrics (no conduction in the dielectrics). Summing the potential
// differences between successive points of the discretisation (electrode
// faces, cell centres, interface faces) from the powered to the ground
// electrode, with D = eps*E on each segment, gives
//
//     j = (d(Vp - Vg)/dt + sum Jc_f L_f/eps_f)/(sum L_f/eps_f)
//
// over the segments f of length L_f: cell centre to cell centre across an
// internal face, cell centre to face at a boundary. Jc_f is the conduction
// current density of the species equations through that face (zero in the
// dielectrics). This holds exactly for the discrete equations when the step
// conserves charge; R_gas measures how far it does not. No difference of E
// is taken.
if (totalCurrent1D)
{
    const scalar dtStep = runTime.deltaT().value();

    const scalarField& magSfG = mesh.magSf().boundaryField()[groundPatchI];
    const scalar area = gSum(magSfG);

    // Unit vector along the axis, towards the ground electrode
    const vector axis = gSum(mesh.Sf().boundaryField()[groundPatchI])/area;

    const scalar Vg = gSum(Phi.boundaryField()[groundPatchI]*magSfG)/area;

    const scalarField& magSfP =
        dielectricRegions[poweredRegionI].magSf().boundaryField()[poweredPatchI];

    const scalar Vp =
        gSum(PhiD[poweredRegionI].boundaryField()[poweredPatchI]*magSfP)
       /gSum(magSfP);

    const scalar dV = Vp - Vg;

    // Charge through each gas face [A], along Sf
    const surfaceScalarField Ic(mspm().netChargeFaceFlux());

    scalar sumLInvEps = 0;
    scalar sumJcLInvEps = 0;

    // Gas: internal faces (permittivity interpolated as in the Laplacian)
    {
        const surfaceScalarField epsf(linearInterpolate(eps));
        const vectorField& C = mesh.C().internalField();
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();
        const vectorField& Sf = mesh.Sf().internalField();
        const scalarField& magSf = mesh.magSf().internalField();

        forAll(own, facei)
        {
            const scalar L = mag((C[nei[facei]] - C[own[facei]]) & axis);
            const scalar jc =
                Ic.internalField()[facei]/magSf[facei]
               *((Sf[facei]/magSf[facei]) & axis);

            sumLInvEps += L/epsf.internalField()[facei];
            sumJcLInvEps += jc*L/epsf.internalField()[facei];
        }

        // Boundary faces: half cells, permittivity of the cell
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];

            if (!p.size())
            {
                continue;
            }

            const unallocLabelList& fc = p.faceCells();
            const vectorField nf(p.nf());

            forAll(fc, k)
            {
                const scalar L = mag((p.Cf()[k] - C[fc[k]]) & axis);
                const scalar jc =
                    Ic.boundaryField()[patchi][k]/p.magSf()[k]
                   *(nf[k] & axis);
                const scalar e = eps.internalField()[fc[k]];

                sumLInvEps += L/e;
                sumJcLInvEps += jc*L/e;
            }
        }
    }

    // Dielectrics: no conduction
    forAll(dielectricRegions, i)
    {
        const fvMesh& dm = dielectricRegions[i];
        const surfaceScalarField epsf(linearInterpolate(epsD[i]));
        const vectorField& C = dm.C().internalField();
        const unallocLabelList& own = dm.owner();
        const unallocLabelList& nei = dm.neighbour();

        forAll(own, facei)
        {
            sumLInvEps +=
                mag((C[nei[facei]] - C[own[facei]]) & axis)
               /epsf.internalField()[facei];
        }

        forAll(dm.boundary(), patchi)
        {
            const fvPatch& p = dm.boundary()[patchi];

            if (!p.size())
            {
                continue;
            }

            const unallocLabelList& fc = p.faceCells();

            forAll(fc, k)
            {
                sumLInvEps +=
                    mag((p.Cf()[k] - C[fc[k]]) & axis)
                   /epsD[i].internalField()[fc[k]];
            }
        }
    }

    // Charge balance of the gas over the step
    const scalar gasCharge = mspm().totalCharge();

    scalar chargeOut = 0;

    forAll(Ic.boundaryField(), patchi)
    {
        chargeOut += sum(Ic.boundaryField()[patchi]);
    }

    const scalar IcondGround = -sum(Ic.boundaryField()[groundPatchI]);

    if (electrodeDVOldValid)
    {
        const scalar j =
            ((dV - electrodeDVOld)/dtStep + sumJcLInvEps)/sumLInvEps;

        Jtot == dimensionedVector("j", Jtot.dimensions(), j*axis);

        forAll(dielectricRegions, i)
        {
            JtotD[i] == dimensionedVector("j", JtotD[i].dimensions(), j*axis);
        }

        // Current from the ground electrode into the gas
        const scalar ItotGround = -j*area;

        const scalar Rgas = (gasCharge - gasChargeOld)/dtStep + chargeOut;

        if (chargeBalanceFile.valid())
        {
            chargeBalanceFile()
                << runTime.value() << tab << dtStep << tab
                << ItotGround << tab << IcondGround << tab
                << ItotGround - IcondGround << tab
                << Rgas << tab << mspm().postSolveCharge()/dtStep << endl;
        }
    }
    else
    {
        Info<< "Total current: effective length g + d/eps_r = "
            << plasmaConstants::epsilon0.value()*sumLInvEps << " m" << endl;
    }

    electrodeDVOld = dV;
    gasChargeOld = gasCharge;
    electrodeDVOldValid = true;
}
