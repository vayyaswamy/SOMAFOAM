 {
	forAll(Phi.boundaryField(), patch)
	{
		//Jnet = mspm().netChargeFlux();
   
   		//Info << Phi.boundaryField()[patch].type() << endl;
     if (Phi.boundaryField()[patch].type() == "dielectricSideWall")
     {
       surfC.boundaryField()[patch] = surfC.oldTime().boundaryField()[patch] + (Jnet.boundaryField()[patch] & mesh.Sf().boundaryField()[patch]) / mesh.magSf().boundaryField()[patch] *runTime.deltaT().value();
     }

		//Info << "Jnet = " << (Jnet.boundaryField()[patchID] & mesh.Sf().boundaryField()[patchID]) / mesh.magSf().boundaryField()[patchID]   << endl;

		//Info << "surfC = " << surfC.boundaryField()[patchID] << endl;

		//Info << "JtotD = " << JtotD[i].boundaryField()[patchID] << endl;
	}
}
