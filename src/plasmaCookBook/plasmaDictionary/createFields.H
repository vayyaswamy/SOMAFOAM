    IOdictionary specieList
    (
        IOobject
        (
            "reactions",    // dictionary name
            runTime.constant(),     
            mesh,                   
            IOobject::MUST_READ,    
            IOobject::NO_WRITE     
        )
    );
    wordList speciesNamesCyp
    (
        specieList.lookup("species")
    );

    IOdictionary activeSpecieList
    (
        IOobject
        (
            "plasmaProperties",    // dictionary name
            runTime.constant(),     
            mesh,                   
            IOobject::MUST_READ,    
            IOobject::NO_WRITE     
        )
    );

	scalar activeSpecies(activeSpecieList.lookupOrDefault<scalar>("activeSpecies",0));

    forAll(speciesNamesCyp, s)
    {
        word currentSpecie = speciesNamesCyp[s];

		volScalarField currentSpecieN
		(
			IOobject
			(
				currentSpecie,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		);
		currentSpecieN.write();

		if(s<activeSpecies)
		{
			if(!Uref)
			{
				volVectorField sFlux
				(
					IOobject
					(
						"F_" + currentSpecie,
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector(0.0, 0.0, 0.0)),
					zeroGradientFvPatchScalarField::typeName
				);
				sFlux.write();
			}
			else
			{
				volVectorField sFlux
				(
					IOobject
					(
						"U_" + currentSpecie,
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector(0.0, 0.0, 0.0)),
					zeroGradientFvPatchScalarField::typeName
				);
				sFlux.write();
			}
		}
    }



	// ------------------------------------------------------------------------- //

	volScalarField Phi
	(
		IOobject
		(
		    "Phi",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0)
	);

	// ------------------------------------------------------------------------- //

	volScalarField p
	(
		IOobject
		(
		    "p",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField T
	(
		IOobject
		(
		    "T",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField Te
	(
		IOobject
		(
		    "Te",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField Tion
	(
		IOobject
		(
		    "Tion",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

    volVectorField eeFlux
    (
        IOobject
        (
            "eeFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
    	dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector(0.0, 0.0, 0.0)),
		zeroGradientFvPatchScalarField::typeName
    );

	Phi.write();

	T.write();

	Te.write();

	Tion.write();

	p.write();

	eeFlux.write();
