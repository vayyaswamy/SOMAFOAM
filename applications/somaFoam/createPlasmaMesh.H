    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    // A dynamic mesh (e.g. adaptive refinement) is used when the case has a
    // constant/dynamicMeshDict; otherwise the mesh is static
    autoPtr<dynamicFvMesh> meshPtr;

    {
        IOobject meshIO
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        );

        IOobject dynamicMeshDictHeader
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::NO_READ
        );

        if (dynamicMeshDictHeader.headerOk())
        {
            meshPtr = dynamicFvMesh::New(meshIO);
        }
        else
        {
            meshPtr.set(new staticFvMesh(meshIO));
        }
    }

    dynamicFvMesh& mesh = meshPtr();
