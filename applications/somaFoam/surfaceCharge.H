{
    forAll(dielectricRegions, i)
    {

        label patchID = mesh.boundaryMesh().findPatchID(rp.dielectricRegionNames()[i]);

        surfC.boundaryField()[patchID] = surfC.oldTime().boundaryField()[patchID] + (Jnet.boundaryField()[patchID] & mesh.Sf().boundaryField()[patchID]) / mesh.magSf().boundaryField()[patchID] *runTime.deltaT().value();


    }
}
