    PtrList<volScalarField> PhiD(dielectricRegions.size());
    PtrList<volScalarField> PhiDOld(dielectricRegions.size());
    PtrList<volVectorField> ED(dielectricRegions.size());
    PtrList<volVectorField> JtotD(dielectricRegions.size());
    PtrList<volScalarField> epsD(dielectricRegions.size());
    PtrList<volVectorField> ddtED(dielectricRegions.size());

    forAll(dielectricRegions, i)
    {
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
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
			    dimensionedVector("E", dimensionSet(0, 0, 0, 1, 0), vector(0.0, 0.0, 0.0))
            )
        );

        JtotD.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "JtotD",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
                dimensionedVector("JtotD", dimensionSet(0, 0, 0, 1, 0), vector(0.0, 0.0, 0.0))
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
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i]
            )
        );
	}

	
    

