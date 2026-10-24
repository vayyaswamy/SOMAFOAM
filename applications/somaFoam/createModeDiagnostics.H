// Per-step quantities that classify the discharge mode (postStepModeDiagnostics.H).
//
// Optional (controlDict modeDiagnostics, default off), 1D with the total
// current of postStepTotalCurrent.H (totalCurrent1D). Written to
// modeDiagnostics/<start time>/modeDiagnostics.dat.
const Switch modeDiagnostics
(
    runTime.controlDict().lookupOrDefault<Switch>("modeDiagnostics", false)
);

const bool modeDiagnostics1D = modeDiagnostics && totalCurrent1D;

if (modeDiagnostics && !totalCurrent1D)
{
    WarningIn(args.executable())
        << "modeDiagnostics needs a 1D plasma mesh with two electrodes"
        << " (see the total current); not written" << endl;
}

label modeDiagElectronI = -1;
label modeDiagDielectricPatchI = -1;
autoPtr<OFstream> modeDiagFile;

if (modeDiagnostics1D)
{
    modeDiagElectronI = mspm().species()["electron"];

    forAll(mesh.boundaryMesh(), patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "regionCouple")
        {
            modeDiagDielectricPatchI = patchi;
        }
    }

    if (modeDiagDielectricPatchI < 0)
    {
        FatalErrorIn(args.executable())
            << "modeDiagnostics: no plasma-dielectric interface patch"
            << exit(FatalError);
    }

    if (Pstream::master())
    {
        const fileName dir(runTime.path()/"modeDiagnostics"/runTime.timeName());

        mkDir(dir);

        modeDiagFile.reset(new OFstream(dir/"modeDiagnostics.dat"));

        modeDiagFile().precision(8);

        modeDiagFile()
            << "# End of each step. y: distance from the dielectric surface"
            << " (y = 0) towards the ground electrode (y = d), as y/d." << nl
            << "# alpha_eff = nu_net/(mu_e |E|), nu_net = d(electron source)/dNe"
            << " (ionization - attachment - recombination frequency)." << nl
            << "# lnM: largest integral of alpha_eff dy over any interval of"
            << " the gap (log of the largest avalanche gain; Raether: ~18-20)."
            << nl
            << "# E_mean = |V_ground - V_dielectric surface|/d;"
            << " E_ground, E_dielectric: field normal to each wall, positive"
            << " towards the gas." << nl
            << "# sigma_+ = e * integral of the positive-ion density dy,"
            << " sigma_e = e * integral of Ne dy [C/m2]." << nl
            << "# t [s]\tdeltaT [s]\tNe_max [1/m3]\ty(Ne_max)/d"
            << "\tNion+_max [1/m3]\tE_max [V/m]\ty(E_max)/d\tE_mean [V/m]"
            << "\tE_ground [V/m]\tE_dielectric [V/m]\tlnM\ty(lnM end)/d"
            << "\tsigma_+ [C/m2]\tsigma_e [C/m2]" << endl;
    }
}
