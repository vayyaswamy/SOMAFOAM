	const dimensionedScalar eCharge = 1.60217662E-19;

	const dimensionedScalar bolzC = 1.38064852E-23;

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
		    IOobject::NO_READ,
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

	// ------------------------------------------------------------------------- //

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

	volScalarField surfC
	(
		IOobject
		(
		    "surfC",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0)
	);
