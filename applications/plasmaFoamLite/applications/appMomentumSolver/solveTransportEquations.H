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
// =============== end of source term ==================================================================









//*****************************************************************************************************
//========================== momentum solver ==========================================================
//ionFlux=0*ionFlux;
forAll (speciesList, ii){
	/*Info << transportModel[ii] << endl;
	velocity[ii] = EField*mobility[ii]*chargeNumber[ii]/abs(chargeNumber[ii]);
	ionFlux += nDensity[ii]*velocity[ii];
	*/
	if (transportModel[ii] == "momentum"){
		if (ii==0){
			pTemperature= 1*eTemp;}
		else
			{pTemperature= gasTemperature*dimensionCorrectionTemp;}

		Info << "Solving momentum equation" << endl;//--------------------
			fvVectorMatrix momEqn
			(
				fvm::ddt(nDensity[ii],velocity[ii]) + 
					fvm::div(fvc::interpolate(nDensity[ii])*(fvc::interpolate(velocity[ii]) & mesh.Sf()),velocity[ii],"div(velocity*nDensity,velocity)")  
							//+fvm::Sp(pSource[ii],velocity[ii])
							+fvm::Sp(nDensity[ii]*eCharge/(mobility[ii]+1e-50*dimensionCorrectionMobility)/speciesMass[ii]/dimensionCorrectionMass,velocity[ii])
							== 
								chargeNumber[ii]/abs(chargeNumber[ii])*eCharge*nDensity[ii]*EField/speciesMass[ii]/dimensionCorrectionMass
								- fvc::grad(kBoltz*pTemperature*nDensity[ii])/speciesMass[ii]/dimensionCorrectionMass
			);
		momEqn.relax();
		momEqn.solve();
	}
}
//======================================================================================================














// Particle continuity  ===============================================================================
//=====================================================================================================
forAll (speciesList, ii){
	if (transportModel[ii] == "driftDiffusion"){
		if (chargeNumber[ii] != 0){
			pVel = EField*mobility[ii]*chargeNumber[ii]/abs(chargeNumber[ii]);}
		else {
			pVel = 0*EField*mobility[ii]*chargeNumber[ii];}
		fvScalarMatrix advEqn
			(   
				fvm::ddt(nDensity[ii])+
					fvm::div((fvc::interpolate(pVel) & mesh.Sf()),nDensity[ii],"div(vel,nDensity)") -
					fvm::laplacian(diffusion[ii], nDensity[ii],"laplacian(difC,nDensity)") + fvm::SuSp(-pSource[ii]/nDensity[ii],nDensity[ii])
			);
		advEqn.relax();
		advEqn.solve();
		(nDensity[ii]).max(cosmicRadiation);
		Flux_[ii]=pVel*nDensity[ii] - diffusion[ii]*fvc::grad(nDensity[ii]);
	}	
	else if (transportModel[ii] == "momentum"){
		fvScalarMatrix advEqn
			(   
				fvm::ddt(nDensity[ii])+
					fvm::div((fvc::interpolate(velocity[ii]) & mesh.Sf()),nDensity[ii],"div(vel,nDensity)") 
					+ fvm::SuSp(-pSource[ii]/nDensity[ii],nDensity[ii])
			);
		advEqn.relax();
		advEqn.solve();
		(nDensity[ii]).max(cosmicRadiation);
		Flux_[ii]=velocity[ii]*nDensity[ii];
	}
}





// energy continuity ====================================================================================
//=======================================================================================================
if (transportModel[0] == "driftDiffusion"){
	jouleHeating= (EField & (nDensity[0]*EField*mobility[0]*chargeNumber[0]/abs(chargeNumber[0])-diffusion[0]*fvc::grad(nDensity[0])))*eCharge;
	totEnergyGen=-eSource-jouleHeating;
	fvScalarMatrix advEEqn
		(   
			fvm::ddt(eDensity) +
				fvm::div((fvc::interpolate((5/3)*EField*mobility[0]*chargeNumber[0]/abs(chargeNumber[0])) & mesh.Sf()),eDensity,"div(vel,eDensity)") -
				fvm::laplacian((5/3)*diffusion[0], eDensity,"laplacian(difC,eDensity)") + fvm::SuSp(-totEnergyGen/eDensity, eDensity)
		);
	advEEqn.relax();
	advEEqn.solve();
	eDensity.max(1e-50);
}
else if (transportModel[0]=="momentum"){
	jouleHeating= (EField & (nDensity[0]*velocity[0]-0*diffusion[0]*fvc::grad(nDensity[0])))*eCharge;
	totEnergyGen=-eSource-jouleHeating;
	fvScalarMatrix advEEqn
		(   
			fvm::ddt(eDensity) +
				fvm::div((fvc::interpolate((5/3)*velocity[0]) & mesh.Sf()),eDensity,"div(vel,eDensity)") -
				fvm::laplacian((5/3)*diffusion[0], eDensity,"laplacian(difC,eDensity)") + fvm::SuSp(-totEnergyGen/eDensity, eDensity)
		);
	advEEqn.relax();
	advEEqn.solve();
	eDensity.max(1e-50);
}
//================ end of energy continuity ==================================================================