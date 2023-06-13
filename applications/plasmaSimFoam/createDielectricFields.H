    PtrList<volScalarField> PhiD(dielectricRegions.size());
    PtrList<volVectorField> ED(dielectricRegions.size());
    PtrList<volScalarField> eps(dielectricRegions.size());
    PtrList<IOdictionary> dielectricProperties(dielectricRegions.size());

    forAll(dielectricRegions, i)
    {
        dielectricProperties.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "dielectricProperties",
                    runTime.constant(),
                    dielectricRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        PhiD.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Phi",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i]
            )
        );

        ED.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "E",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
			    dimensionedVector("E", dimensionSet(0, 0, 0, 1, 0), vector(0.0, 0.0, 0.0))
            )
        );

        eps.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "eps",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
			    dimensionedScalar("eps", dimensionSet(0, 0, 0, 1, 0), dielectricProperties[i].lookup("epsilon"))
            )
        );
	}

	const fvPatchList& patcha = mesh.boundary();

	volScalarField Ft
	(
		IOobject
		(
		    "Ft",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
	   dimensionedScalar("Ft", dimLength/dimTime, 0.0)
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
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

