// Adaptive mesh refinement, part 1 (before the mesh update): data for the
// transfer of the densities and the electron temperature to refined cells.
//
// The mesh update gives all the cells that a refined cell is split into
// the value of that cell. For the slowly evolving densities that staircase
// persists and disturbs the field. Here the centres of the cells and the
// gradients are stored in fields that the mesh update maps in the same way,
// so that afterwards every new cell knows the centre and the gradient of
// the cell it came from (see prolongApply.H).

autoPtr<volVectorField> prolongCentres;
PtrList<volVectorField> prolongGradLogN(composition.species().size());
autoPtr<volVectorField> prolongGradTe;

if (!isA<staticFvMesh>(mesh))
{
    prolongCentres.reset
    (
        new volVectorField
        (
            IOobject
            (
                "prolongCentres",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimLength, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );
    prolongCentres().internalField() = mesh.C().internalField();

    forAll(prolongGradLogN, i)
    {
        const volScalarField& Ni = mspm().N(i);

        prolongGradLogN.set
        (
            i,
            new volVectorField
            (
                "prolongGradLogN_" + composition.species()[i],
                fvc::grad(Ni)/(Ni + dimensionedScalar("small", Ni.dimensions(), SMALL))
            )
        );
    }

    prolongGradTe.reset
    (
        new volVectorField("prolongGradTe", fvc::grad(thermo.Te()))
    );
}
