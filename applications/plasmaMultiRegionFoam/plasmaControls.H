	IOdictionary plasmaProperties
	(
		IOobject
		(
		    "plasmaProperties",
		    runTime.constant(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	);

	IOdictionary electroMagnetics
	(
		IOobject
		(
		    "electroMagnetics",
		    runTime.constant(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	);

	const word poissonEquationSolver(electroMagnetics.lookup("poissonEquationSolver"));

	if ((poissonEquationSolver != "semiImplicit") && (poissonEquationSolver != "explicit") && (poissonEquationSolver != "linearized"))
	{
		FatalError << "in electroMagnetics, unknown poissonEquationSolver type " << nl
		    << poissonEquationSolver << abort(FatalError);
	}

	const word emcModelRef(electroMagnetics.lookup("emcModel"));

	Switch restartCapabale(runTime.controlDict().lookup("restartCapable"));

    const dictionary& plasmaControl =
        mesh.solutionDict().subDict("plasma");

    int nCorr(readInt(plasmaControl.lookup("nCorrectors")));

    scalar convergenceTolerance(readScalar(plasmaControl.lookup("poissonTolerance")));
