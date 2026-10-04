    // decomposePar does not write the values of the cells on the other
    // processor on the processor patches: fetch them before they are used,
    // so that the first time steps of a parallel run match a serial run
    if (Pstream::parRun())
    {
        HashTable<const volScalarField*> sFields
        (
            mesh.objectRegistry::lookupClass<volScalarField>()
        );

        forAllConstIter(HashTable<const volScalarField*>, sFields, iter)
        {
            const_cast<volScalarField&>(*iter()).boundaryField()
                .evaluateCoupled();
        }

        HashTable<const volVectorField*> vFields
        (
            mesh.objectRegistry::lookupClass<volVectorField>()
        );

        forAllConstIter(HashTable<const volVectorField*>, vFields, iter)
        {
            const_cast<volVectorField&>(*iter()).boundaryField()
                .evaluateCoupled();
        }
    }
