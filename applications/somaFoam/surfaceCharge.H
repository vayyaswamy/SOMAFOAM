{
	forAll(dielectricRegions, i)
	{
		//Jnet = mspm().netChargeFlux();

		label patchID = mesh.boundaryMesh().findPatchID(rp.dielectricRegionNames()[i]);			

		surfC.boundaryField()[patchID] = surfC.oldTime().boundaryField()[patchID] + (Jnet.boundaryField()[patchID] & mesh.Sf().boundaryField()[patchID]) / mesh.magSf().boundaryField()[patchID] *runTime.deltaT().value();

		//surfC.boundaryField()[patchID] = 1e12*8.85e-12;

		//Info << "Jnet = " << (Jnet.boundaryField()[patchID] & mesh.Sf().boundaryField()[patchID]) / mesh.magSf().boundaryField()[patchID]   << endl;

		//Info << "surfC = " << surfC.boundaryField()[patchID] << endl;

		//Info << "JtotD = " << JtotD[i].boundaryField()[patchID] << endl;
	}
}
