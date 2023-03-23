	const dimensionedScalar eCharge = 1.60217662E-19;
	const dimensionedScalar bolzC = 1.38064852E-23;

	// dictionary lookup in sys/controlDic for courant number parameters 
	scalar MaxCo(readScalar(runTime.controlDict().lookup("MaxCo")));
	scalar deltaTMax(readScalar(runTime.controlDict().lookup("deltaTMax")));
	scalar deltaTMin(readScalar(runTime.controlDict().lookup("deltaTMin")));

	// controls for plasma simulation and species utilized 
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
	const word emcModelRef(electroMagnetics.lookup("emcModel"));

	// dictionary lookup in const/electroMagnetics for spatio-temporal controls
	const word refSpatioTemporal(electroMagnetics.lookup("spatioTemporal"));

	// determine if the courant criterion to be utilized or not
	const word courantNumber(runTime.controlDict().lookup("courantNumber"));

	if ((solutionDomain != "plasmaDielectric") && (solutionDomain != "plasma"))
	{
		FatalError << "in electroMagnetics, unknown solutionDomain type " << nl
		    << solutionDomain << abort(FatalError);
	}

	// error checking for model solving the Poisson's equation
	if ((poissonEquationSolver != "semiImplicit") && (poissonEquationSolver != "explicit") && (poissonEquationSolver != "linearized"))
	{
		FatalError << "in electroMagnetics, unknown poissonEquationSolver type " << nl
		    << poissonEquationSolver << abort(FatalError);
	}

	// dictionary lookup in sys/controlDict for user-defined simulation control  
	Switch restartCapabale(runTime.controlDict().lookup("restartCapable"));

	// initialize chemitry modules and reference different components
	autoPtr<psiChemistryModel> pChemistry 
	(
		psiChemistryModel::New(mesh) 
	);

	psiChemistryModel& chemistry = pChemistry();
	hsCombustionThermo& thermo = chemistry.thermo();
	basicMultiComponentMixture& composition = thermo.composition();
	PtrList<volScalarField>& Y = composition.Y();
	PtrList<scalar> W_(composition.species().size());

	// ------------------------------------------------------------------------- //

	// mesh scalar of voltage in the plasma
	volScalarField Phi
	(
		IOobject
		(
		    "Phi",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		),
		mesh
	);

	mesh.schemesDict().setFluxRequired(Phi.name());

	volScalarField epsilon0
	(
		IOobject
		(
		    "epsilon",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("epsilon", dimensionSet(0, 0, 0, 1, 0), 8.85418782E-12)
	);

	volScalarField epsiloni
	(
		IOobject
		(
		    "epsiloni",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("epsiloni", dimensionSet(0, 0, 0, 1, 0), 0)
	);

	volVectorField E
	(
		IOobject
		(
		    "E",
		    runTime.timeName(),
		    mesh,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

	volVectorField gradTe
	(
		IOobject
		(
		    "gradTe",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

	// electron conductivity within the system
	volScalarField kappa_e
	(
		IOobject
		(
		    "kappa_electron",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0)
	);

	volVectorField ddtE
	(
		IOobject
		(
		    "ddtE",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

	// for electromagnetic modules with the plasmaExtendedPotential BC, this is required 
	IOdictionary voltageDict
	(
		IOobject
		(
		    "voltageDict",
		    runTime.constant(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		)
	);

	volScalarField rho
	(
		IOobject
		(
		    "rho",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		thermo.rho()
	);

	volVectorField ionFlux
	(
		IOobject
		(
		    "ionFlux",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

  	volVectorField conductionCurrent
    (
    	IOobject
    	(
    		"conductionCurrent",
    		runTime.timeName(),
    		mesh, 
    		IOobject::NO_READ,
    		IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("conductionCurrent", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

  	volVectorField totalCurrent
    (
    	IOobject
    	(
    		"totalCurrent",
    		runTime.timeName(),
    		mesh, 
    		IOobject::NO_READ,
    		IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("totalCurrent", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

    // volumetric scalar field for sources and sinks of electrons in the system. 
    // Calculation obtained by obtaining the divergence of electron flux
	volScalarField electronDiv
	(
		IOobject
		(
			"electronDiv",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0)
	);

	// hastables for species in the system 
	multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

	forAll(Y, i)
	{
		fields.add(Y[i]);

		mesh.schemesDict().setFluxRequired(Y[i].name());

		W_.set
		(
		    i,
		    new scalar
		    (
		        dynamic_cast<const reactingMixture<constGasThermoPhysics>&>
		            (thermo).speciesData()[i].W()
		    )
		);
	}

	// plasma module defined 
	autoPtr<multiSpeciesPlasmaModel> mspm
	(
		multiSpeciesPlasmaModel::New
		(
		    thermo
		)
	);

	// initialized plasma module
	mspm().input();

	// initializes energy module
	plasmaEnergyModel pem(thermo, mspm(), E);

	// initializes electromagnetic module 
    autoPtr<emcModel> emc
	(
		emcModel::New
		(
		    electroMagnetics,
		    mspm(),
		    E,
			runTime
		)
    );

	// initializes spatioTemporal Module 
	autoPtr<spatioTemporal> spatioTemp
	(
		spatioTemporal::New
		(
			electroMagnetics,
		    electroMagnetics,
		    mspm(),
		    E,
			runTime,
			mesh
		)
    );

