    regionProperties rp(runTime);

    PtrList<fvMesh> dielectricRegions(rp.dielectricRegionNames().size());

    forAll(rp.dielectricRegionNames(), i)
    {
        dielectricRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp.dielectricRegionNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
