//***************************************************************************************************************************
// CALCULATING THE NET CHARGE ALL ALL GRID POINTS INSIDE THE SIMULATION DOMAIN **********************************************
//===========================================================================================================================
netCharge=0*netCharge;
forAll(speciesList, ii){
		netCharge= netCharge + eCharge*chargeNumber[ii]*nDensity[ii];
}
//***************************************************************************************************************************
//do while loop to test the convergence *************************************************************************************
//===========================================================================================================================
float maxValue=0.0;
int nIterations1=0;
do{
		EPotential.storePrevIter();

		//********************************************************************************************************************
		// Dirichlet Boundary condition (applied at plasma side bounraries) **************************************************
		//====================================================================================================================
		forAll (dielectricList, dl){
			//Info << "Dirichlet ma aayo" << endl;
			label patchID = mesh.boundaryMesh().findPatchID(dielectricList[dl]);
			const polyPatch plasma_patch= mesh.boundaryMesh()[patchID];
			const fvPatchField<scalar> plasma_patch_field = EPotential.boundaryField()[patchID];
			forAll (dielectricList, jk){
				label patchID1 = (dielectricMesh[jk]).boundaryMesh().findPatchID(dielectricList[dl]);
				if (patchID1>-1){
					const polyPatch dielectric_patch= dielectricMesh[jk].boundaryMesh()[patchID1];
					const fvPatchField<scalar> dielectric_patch_field = EPotentialD[jk].boundaryField()[patchID1];
					patchToPatchInterpolation interpolator(plasma_patch,dielectric_patch); 	
					scalarField plasmaToDielectric = interpolator.faceInterpolate(plasma_patch_field);
					forAll (dielectricMesh[jk].boundary()[patchID1],faceII){
						//EPotentialD[jk].boundaryField()[patchID1][faceII]=EPotentialD[jk].boundaryField()[patchID1][faceII]+0.2*(plasmaToDielectric[faceII]-EPotentialD[jk].boundaryField()[patchID1][faceII]);	
						 EPotentialD[jk].boundaryField()[patchID1][faceII]=plasmaToDielectric[faceII];	
					}
				}
			}
		}

		//********************************************************************************************************************
		// Solving the poisson equation for the plasma region ****************************************************************
		//====================================================================================================================
		fvScalarMatrix poissonPlasmaEquation
						(
							fvm::laplacian(EPotential) + (netCharge) /epsilon //explicit term
						);
		forAll(speciesList, kp){ //adding from future timestep (semi-implicit)
			poissonPlasmaEquation += 
					chargeNumber[kp]*eCharge*runTime.deltaT()/epsilon*fvm::laplacian(nDensity[kp]*mobility[kp], EPotential,"laplacian(nmu,E)")
					+ chargeNumber[kp]*eCharge*runTime.deltaT()/epsilon*fvc::laplacian(diffusion[kp],nDensity[kp],"laplacian(difC,nDensity)");
		}
		poissonPlasmaEquation.relax();
		poissonPlasmaEquation.solve();
		EField=-fvc::grad(EPotential);
		EField.correctBoundaryConditions();



		//********************************************************************************************************************
		// Neumann Boundary condition (applied at dielectric side boundaries) ************************************************
		//====================================================================================================================
		forAll (dielectricList, dl){
			//Info << "Neumann ma aayo" << endl;
			label patchID = mesh.boundaryMesh().findPatchID(dielectricList[dl]);
			const polyPatch plasma_patch= mesh.boundaryMesh()[patchID];
			const fvPatchField<scalar> plasma_patch_field = EPotential.boundaryField()[patchID];
			forAll (dielectricList, jk){
				label patchID1 = (dielectricMesh[jk]).boundaryMesh().findPatchID(dielectricList[dl]);
				if (patchID1>-1){
					const polyPatch dielectric_patch= dielectricMesh[jk].boundaryMesh()[patchID1];
					vectorField dielectric_patch_field = 0*EFieldD[jk].boundaryField()[patchID1];
					vectorField normmD= dielectricMesh[jk].boundary()[patchID1].nf();
					forAll (dielectricMesh[jk].boundary()[patchID1],faceDi){
						dielectric_patch_field[faceDi]=(EPotentialD[jk][dielectricMesh[jk].boundary()[patchID1].faceCells()[faceDi]]-
															EPotentialD[jk].boundaryField()[patchID1][faceDi])*
															(dielectricMesh[jk].surfaceInterpolation::deltaCoeffs()).boundaryField()[patchID1][faceDi]
															*normmD[faceDi];
					 } 
					const fvPatchField<scalar> dielectric_patch_field_epsilon = epsilonD[jk].boundaryField()[patchID1];
					patchToPatchInterpolation interpolator(dielectric_patch,plasma_patch); 	
					vectorField dielectricToPlasma = interpolator.faceInterpolate(dielectric_patch_field);

					scalarField dielectricToPlasmaEpsilon = interpolator.faceInterpolate(dielectric_patch_field_epsilon);					
					vectorField normm= mesh.boundary()[patchID].nf();
					forAll (mesh.boundary()[patchID],faceII){
							EPotential.boundaryField()[patchID][faceII]=EPotential[mesh.boundary()[patchID].faceCells()[faceII]]
							-(  (dielectricToPlasma[faceII] & normm[faceII])*dielectricToPlasmaEpsilon[faceII]
							+0*surfC[mesh.boundary()[patchID].faceCells()[faceII]])
							/((mesh.surfaceInterpolation::deltaCoeffs()).boundaryField()[patchID][faceII]*epsilon[mesh.boundary()[patchID].faceCells()[faceII]]);
					}
				}
			}
		}

		//********************************************************************************************************************
		// Solving the poisson equation for all dielectric region ************************************************************
		//====================================================================================================================
		forAll (dielectricList,ii){
			EPotentialD[ii].storePrevIter();
			fvScalarMatrix poissonDielectricEquation
				(
					fvm::laplacian(EPotentialD[ii], "laplacian(EPotential)")
				);
				poissonDielectricEquation.relax();
				poissonDielectricEquation.solve();
				EFieldD[ii]=-fvc::grad(EPotentialD[ii]);
				EFieldD[ii].correctBoundaryConditions();
				//Info << "Dielectric lai solve gardai cha" << endl;
		}


		//*********************************************************************************************************************
		// stability criteria check. ******************************************************************************************
		//=====================================================================================================================
		if (dielectricList.size()>0){
			//Info << "convergence check hudai" << endl;
			maxValue=gMax(mag(EPotential.internalField()-EPotential.prevIter().internalField()));
			forAll (dielectricList,dlec){
				float maxValue1 = gMax(mag(EPotentialD[dlec].internalField()-EPotentialD[dlec].prevIter().internalField()));
				if (maxValue1>maxValue)
					maxValue=maxValue1;	
			}//====================================================================================================================
			//Info << "MaxValue::" << maxValue << "convergenceDelta::" << convergenceDelta << endl;
		} 
		nIterations1++;
}
while (maxValue>convergenceDelta);
Info << "Solved Poisson ======================================== iterations"<< nIterations1 <<endl;