forAll (speciesList,ik){		
	nDensity[ik].correctBoundaryConditions();
	eDensity.correctBoundaryConditions();
	double exitProbability=0;
	// =================================== FOR METALLIC WALLS ====================================================	
	forAll(metallicBoundary, ii)
	{
		label patchID = mesh.boundaryMesh().findPatchID(metallicBoundary[ii]);
		forAll (mesh.boundary()[patchID],faceI){
				particleFlux=0*particleFlux;
				energyFlux=0*energyFlux;
				vectorField normm= mesh.boundary()[patchID].nf();
				exitProbability=1;//Foam::exp(-1.6e-19*abs(EPotential[mesh.boundary()[patchID].faceCells()[faceI]]-EPotential.boundaryField()[patchID][faceI])/
						//		(1.38e-23*eTemp[mesh.boundary()[patchID].faceCells()[faceI]]));
				mobility[ii].boundaryFieldRef()[patchID][faceI]=0;
				// thermal flux ---------------------------------
				if (ik==0){
					particleFlux = 0*1/2*Foam::sqrt((8*1.38e-23*eTemp[mesh.boundary()[patchID].faceCells()[faceI]]/((22/7)*9.1e-31)))*
									nDensity[ik][mesh.boundary()[patchID].faceCells()[faceI]]*normm[faceI]*exitProbability;
					energyFlux = 0*5/6*Foam::sqrt((8*1.38e-23*eTemp[mesh.boundary()[patchID].faceCells()[faceI]]/((22/7)*9.1e-31)))*
									eDensity[mesh.boundary()[patchID].faceCells()[faceI]]*normm[faceI]*exitProbability;
				}

				if ((normm[faceI] & ((velocity[ik])[mesh.boundary()[patchID].faceCells()[faceI]])) > 0.0) 
					{// moving towards boundary.
						particleFlux += nDensity[ik][mesh.boundary()[patchID].faceCells()[faceI]]*velocity[ik][mesh.boundary()[patchID].faceCells()[faceI]];
						mobility[ik].boundaryFieldRef()[patchID][faceI] = mobility[ik][mesh.boundary()[patchID].faceCells()[faceI]];
						if (ik==0)
							energyFlux += eDensity[mesh.boundary()[patchID].faceCells()[faceI]]*eVelocity[mesh.boundary()[patchID].faceCells()[faceI]];
					}
                                pFlux[ik].boundaryFieldRef()[patchID][faceI] = particleFlux/nDensity[ik][mesh.boundary()[patchID].faceCells()[faceI]];
				eFlux.boundaryFieldRef()[patchID][faceI]=energyFlux/eDensity[mesh.boundary()[patchID].faceCells()[faceI]];
		}
	}
}// ======================================================================================================================

	
