PtrList<volScalarField> EPotentialD(dielectricList.size());
PtrList<volVectorField> EFieldD(dielectricRegions.size());
PtrList<volVectorField> EFieldDManual(dielectricRegions.size());
PtrList<volScalarField> epsilonD(dielectricList.size());
PtrList<volScalarField> surfCD(dielectricList.size());
forAll(dielectricList, i)// only when there are dielectric layers
{
    Info << i << dielectricList[i] <<endl ;
    EPotentialD.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "EPotential",
                runTime.timeName(),
                dielectricMesh[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricMesh[i]
        )
    );
    

    surfCD.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "surfC",
                runTime.timeName(),
                dielectricMesh[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricMesh[i]
        )
    );


    EFieldD.set
    (
        i,
        new volVectorField 
        (
            IOobject
            (
                "EField",
                runTime.timeName(),
                dielectricMesh[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricMesh[i],
            dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );

       EFieldDManual.set
    (
        i,
        new volVectorField 
        (
            IOobject
            (
                "EFieldManual",
                runTime.timeName(),
                dielectricMesh[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dielectricMesh[i],
            dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)//,
            //zeroGradientFvPatchVectorField::typeName
        )
    );

    epsilonD.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime.timeName(),
                dielectricMesh[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricMesh[i]
        )
    );

}