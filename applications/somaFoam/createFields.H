	const dimensionedScalar eCharge = 1.60217662E-19;

	const dimensionedScalar bolzC = 1.38064852E-23;

	scalar MaxCo(readScalar(runTime.controlDict().lookup("MaxCo")));

	scalar deltaTMax(readScalar(runTime.controlDict().lookup("deltaTMax")));

	scalar deltaTMin(readScalar(runTime.controlDict().lookup("deltaTMin")));

		

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

	const word solutionDomain(electroMagnetics.lookup("solutionDomain"));

	const word poissonEquationSolver(electroMagnetics.lookup("poissonEquationSolver"));

	const word emcModelRef(electroMagnetics.lookup("emcModel"));

	if ((solutionDomain != "plasmaDielectric") && (solutionDomain != "plasma"))
	{
		FatalError << "in electroMagnetics, unknown solutionDomain type " << nl
		    << solutionDomain << abort(FatalError);
	}

	if ((poissonEquationSolver != "semiImplicit") && (poissonEquationSolver != "explicit") && (poissonEquationSolver != "linearized"))
	{
		FatalError << "in electroMagnetics, unknown poissonEquationSolver type " << nl
		    << poissonEquationSolver << abort(FatalError);
	}

	Switch restartCapabale(runTime.controlDict().lookup("restartCapable"));

	Info << "Chemistry step ? " << endl;

	autoPtr<psiChemistryModel> pChemistry 
	(
		psiChemistryModel::New(mesh) 
	);

	//Info << "Chemistry steo done ? " << endl;

	//Info << "Step 1" << endl;

	psiChemistryModel& chemistry = pChemistry();

	//Info << "Step 2" << endl;

	hsCombustionThermo& thermo = chemistry.thermo();

	//Info << "Step 3" << endl;

	basicMultiComponentMixture& composition = thermo.composition();

	//Info << "Step 4" << endl;

	PtrList<volScalarField>& Y = composition.Y();

	PtrList<scalar> W_(composition.species().size());

	// ------------------------------------------------------------------------- //

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

	volScalarField eps
	(
		IOobject
		(
			"epsilon",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("epsilon", dimensionSet(0, 0, 0, 1, 0), 8.85418782E-12)
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

	volVectorField Jnet
	(
		IOobject
		(
		    "Jnet",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

	volVectorField Jtot
	(
		IOobject
		(
		    "Jtot",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
	);

	

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

	autoPtr<multiSpeciesPlasmaModel> mspm
	(
		multiSpeciesPlasmaModel::New
		(
		    thermo
		)
	);

	mspm().input();

	plasmaEnergyModel pem(thermo, mspm(), E);

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

