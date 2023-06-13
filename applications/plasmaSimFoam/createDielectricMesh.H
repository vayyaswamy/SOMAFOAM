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

    forAll(dielectricRegions, i)
    {
		const polyPatchList& patches = dielectricRegions[i].boundaryMesh();

		forAll (patches, patchI)
		{
		    if (isType<regionCouplePolyPatch>(patches[patchI]))
		    {
		        const regionCouplePolyPatch& rcp =
		            refCast<const regionCouplePolyPatch>(patches[patchI]);

		        rcp.attach();
		    }
		}
		dielectricRegions[i].surfaceInterpolation::movePoints();
	}

