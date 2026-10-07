// Run the electronBoltzmann utility (system/boltzmannDict) before the
// transport and rate tables are read, if constant/plasmaProperties has
// "electronBoltzmann yes;". The master runs it; the others wait.
{
    label status = 0;

    if (Pstream::master())
    {
        const fileName casePath(args.rootPath()/args.globalCaseName());

        IFstream is(casePath/"constant"/"plasmaProperties");
        const dictionary properties(is);

        if (properties.lookupOrDefault<Switch>("electronBoltzmann", false))
        {
            Info<< "Running electronBoltzmann (output in "
                << "log.electronBoltzmann)" << nl << endl;

            status = Foam::system
            (
                "electronBoltzmann -case \"" + casePath + "\" > \""
              + casePath/"log.electronBoltzmann" + "\" 2>&1"
            );
        }
    }

    reduce(status, maxOp<label>());

    if (status != 0)
    {
        FatalErrorIn(args.executable())
            << "electronBoltzmann failed: see log.electronBoltzmann"
            << exit(FatalError);
    }
}
