{
	forAll(dielectricRegions, i)
	{
		label patchID = mesh.boundaryMesh().findPatchID(rp.dielectricRegionNames()[i]);			

		Ft.boundaryField()[patchID] = F.boundaryField()[patchID].patchInternalField() & mesh.Sf().boundaryField()[patchID] / mesh.V()[patchID];

		forAll(patcha, patchi)
		{
			const fvPatch& currPatch = patcha[patchID];

			forAll(currPatch, facei)
			{
				label faceCelli = currPatch.faceCells()[facei];
				Ft[faceCelli] = Ft.boundaryField()[patchID][facei];            
			}
		}
	}

    const volScalarField& Fn = Ft;

	solve (fvm::ddt(surfC) - eCharge*Fn);
}
