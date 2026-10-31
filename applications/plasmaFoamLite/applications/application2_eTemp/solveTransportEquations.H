//*****************************************************************************************************
//===== LETS CALCULATE THE VELOCITY OF ALL THE CHARGED SPECIES ========================================
ionFlux=0*ionFlux;
forAll (speciesList, ii){
	if (chargeNumber[ii]!= 0){
		velocity[ii] = EField*mobility[ii]*chargeNumber[ii]/abs(chargeNumber[ii]);
	}
	else{
		velocity[ii] = 0*gasVelocity;
	}
	if (ii != 0)
	ionFlux += nDensity[ii]*velocity[ii];
}

//*****************************************************************************************************
//============ ACCOUNTING FOR THE PRODUCTION OF NEW PARTICLES DUE TO CHEMICAL REACTIONS ===============
eSource=0*eSource;
for (int ii=0; ii<reactionCount ; ii++){
	reactionSource[ii]=reaction[ii]/dimensionCorrectionTime;
	for (int jj=0; jj<numberOfSpecies; jj++){
			reactionSource[ii]=reactionSource[ii]*pow(nDensity[jj],reactantMatrix[ii][jj])
										/pow(dimensionCorrection,reactantMatrix[ii][jj]);
		}
	eSource= eSource + reactionEnergy[ii]*dimensionCorrectionEnergy*reactionSource[ii]/dimensionCorrection;
}

//******************************************************************************************************
// =========================== SOURCE TERM =============================================================
forAll (speciesList, ii){
	pSource[ii]=0*pSource[ii];
	for (int jj=0;jj<reactionCount;jj++){
		pSource[ii]=pSource[ii]+(productMatrix[jj][ii]-reactantMatrix[jj][ii])*reactionSource[jj];
	}
}
// ***************************** JOULE HEATING TERM ***************************************************
//=====================================================================================================
jouleHeating= (EField & (nDensity[0]*velocity[0]-diffusion[0]*fvc::grad(nDensity[0])))*eCharge;
totEnergyGen=-eSource-jouleHeating;

// CONTINUITY =========================================================================================
//=====================================================================================================
forAll (speciesList, ii){
	//nDensity[ii].correctBoundaryConditions();
	fvScalarMatrix advEqn
		(   
			fvm::ddt(nDensity[ii])+
				fvm::div((fvc::interpolate(velocity[ii]) & mesh.Sf()),nDensity[ii],"div(vel,nDensity)") -
				fvm::laplacian(diffusion[ii], nDensity[ii],"laplacian(difC,nDensity)") + fvm::SuSp(-pSource[ii]/nDensity[ii],nDensity[ii])
		);
	advEqn.relax();
	advEqn.solve();
	(nDensity[ii]).max(cosmicRadiation);	
	// energy continuity =============================================
	if (ii==0){
		//eDensity.correctBoundaryConditions();
		fvScalarMatrix advEEqn
			(   
				/*fvm::ddt(eDensity) +
					fvm::div((fvc::interpolate((5/3)*velocity[0]) & mesh.Sf()),eDensity,"div(vel,eDensity)") -
					fvm::laplacian((5/3)*diffusion[0], eDensity,"laplacian(difC,eDensity)") + fvm::SuSp(-totEnergyGen/eDensity, eDensity)
				*/
				fvm::ddt((3/2)*kBoltz*nDensity[0],eTemp)
				+ fvm::div((fvc::interpolate((5/2)*(eVelocity*nDensity[0]*kBoltz)) & mesh.Sf()),eTemp,"div(vel,eDensity)") 
				- fvm::laplacian((5/2)*diffusion[0]*nDensity[0]*kBoltz, eTemp,"laplacian(difC,eDensity)") 
				+ fvm::SuSp(-totEnergyGen/eTemp, eTemp)
			);
		advEEqn.relax();
		advEEqn.solve();
		eDensity.max(300);
	}
}// =================================================================
