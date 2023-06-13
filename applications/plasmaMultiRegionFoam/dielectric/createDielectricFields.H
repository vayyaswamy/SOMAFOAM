    PtrList<volScalarField> PhiD(dielectricRegions.size());

    PtrList<volVectorField> ED(dielectricRegions.size());

    PtrList<volScalarField> epsD(dielectricRegions.size());

    PtrList<volScalarField> sigmaD(dielectricRegions.size());

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
			    dimensionedVector("E", dimensionSet(0, 0, 0, 1, 0), vector::zero)
            )
        );

        epsD.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "epsilon",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
			    dimensionedScalar("epsilon", dimensionSet(0, 0, 0, 1, 0), dielectricProperties[i].lookup("epsilon"))
            )
        );

        sigmaD.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "surfC",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
			    dimensionedScalar("surfC", dimensionSet(0, 0, 0, 1, 0), 0.0)
            )
        );
	}

    // hardcoded IOdictionary for cases where only a input 
    // source is considered
    IOdictionary voltageDict
    (
        IOobject
        (
            "voltageDict",
            runTime.constant(),
            dielectricRegions[0],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

