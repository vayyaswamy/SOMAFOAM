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

	// dictionary lookup in const/electroMagnetics for electrical controls  
	const word solutionDomain(electroMagnetics.lookup("solutionDomain"));

	const word poissonEquationSolver(electroMagnetics.lookup("poissonEquationSolver"));

	if ((poissonEquationSolver != "semiImplicit") && (poissonEquationSolver != "explicit") && (poissonEquationSolver != "linearized"))
	{
		FatalError << "in electroMagnetics, unknown poissonEquationSolver type " << nl
		    << poissonEquationSolver << abort(FatalError);
	}

	const word emcModelRef(electroMagnetics.lookup("emcModel"));

	Switch restartCapabale(runTime.controlDict().lookup("restartCapable"));

	// dictionary lookup in const/electroMagnetics for spatio-temporal controls
	const word refSpatioTemporal(electroMagnetics.lookup("spatioTemporal"));

	// determine if the courant criterion to be utilized or not
	const word courantNumber(runTime.controlDict().lookup("courantNumber"));

    const dictionary& plasmaControl =
        mesh.solutionDict().subDict("plasma");

    int nCorr(readInt(plasmaControl.lookup("nCorrectors")));

    scalar convergenceTolerance(readScalar(plasmaControl.lookup("poissonTolerance")));

	// dictionary lookup in sys/controlDic for courant number parameters 
	scalar MaxCo(readScalar(runTime.controlDict().lookup("MaxCo")));
	scalar deltaTMax(readScalar(runTime.controlDict().lookup("deltaTMax")));
	scalar deltaTMin(readScalar(runTime.controlDict().lookup("deltaTMin")));
