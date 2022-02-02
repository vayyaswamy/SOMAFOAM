{
    const polyPatchList& patches = mesh.boundaryMesh();

    forAll (patches, patchI)
    {
        if (isType<regionCouplePolyPatch>(patches[patchI]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchI]);

            rcp.attach();
        }
    }
    mesh.surfaceInterpolation::movePoints();
}
