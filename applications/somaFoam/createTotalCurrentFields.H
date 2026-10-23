// Total current in 1D (postStepTotalCurrent.H).
//
// Used when the plasma mesh is 1D, with one fixed-potential electrode on the
// plasma mesh and one on a dielectric region; otherwise Jtot keeps the
// expression of plasmaEqn.H.
const Switch consistentTotalCurrent
(
    runTime.controlDict().lookupOrDefault<Switch>("consistentTotalCurrent", true)
);

// Electrode patches: Phi fixed, not coupled, not the plasma-dielectric
// interface
label groundPatchI = -1;
label poweredRegionI = -1;
label poweredPatchI = -1;
bool totalCurrent1D =
    consistentTotalCurrent && mesh.nGeometricD() == 1 && !Pstream::parRun();

if (totalCurrent1D)
{
    label nElectrodes = 0;

    forAll(Phi.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = Phi.boundaryField()[patchi];

        if
        (
            pp.fixesValue() && !pp.coupled()
         && pp.patch().type() != "regionCouple"
         && returnReduce(pp.size(), sumOp<label>()) > 0
        )
        {
            groundPatchI = patchi;
            nElectrodes++;
        }
    }

    forAll(dielectricRegions, i)
    {
        forAll(PhiD[i].boundaryField(), patchi)
        {
            const fvPatchScalarField& pp = PhiD[i].boundaryField()[patchi];

            if
            (
                pp.fixesValue() && !pp.coupled()
             && pp.patch().type() != "regionCouple"
             && returnReduce(pp.size(), sumOp<label>()) > 0
            )
            {
                poweredRegionI = i;
                poweredPatchI = patchi;
                nElectrodes++;
            }
        }
    }

    if (nElectrodes != 2 || groundPatchI < 0 || poweredPatchI < 0)
    {
        totalCurrent1D = false;
    }
}

if (totalCurrent1D)
{
    Info<< "Total current: 1D, from the applied voltage and the conduction"
        << " current (electrodes " << mesh.boundaryMesh()[groundPatchI].name()
        << " and " << dielectricRegions[poweredRegionI].name() << "/"
        << dielectricRegions[poweredRegionI].boundaryMesh()[poweredPatchI].name()
        << ")" << endl;
}
else
{
    Info<< "Total current: conduction + eps*dE/dt between steps" << endl;
}

// Electrode voltage difference and gas charge at the end of the previous
// step
scalar electrodeDVOld = 0;
scalar gasChargeOld = 0;
bool electrodeDVOldValid = false;

// Per-step output: electrode currents and charge balance of the gas
autoPtr<OFstream> chargeBalanceFile;

if (totalCurrent1D && Pstream::master())
{
    const fileName dir(runTime.path()/"chargeBalance"/runTime.timeName());

    mkDir(dir);

    chargeBalanceFile.reset(new OFstream(dir/"chargeBalance.dat"));

    chargeBalanceFile().precision(10);

    chargeBalanceFile()
        << "# End of each step. Currents at the ground electrode, positive"
        << " from the electrode into the gas." << nl
        << "# I_total from the electrode voltages and the face conduction"
        << " fluxes (exact in 1D if charge is conserved)." << nl
        << "# R_gas = dQ_gas/dt + charge flux out through all gas walls:"
        << " zero if the step conserves charge." << nl
        << "# R_postSolve: part of R_gas from relaxation and limits applied"
        << " after the species solution." << nl
        << "# t [s]\tdeltaT [s]\tI_total [A]\tI_conduction [A]"
        << "\tI_displacement [A]\tR_gas [A]\tR_postSolve [A]" << endl;
}
