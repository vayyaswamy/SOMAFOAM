//*****************************************************************************************************
//============================ PARTICLE CONTINUITY  ===================================================
//=====================================================================================================
gradTe= fvc::grad(eTemp);
totalElectronFlux=0*totalElectronFlux;
totalIonFlux=0*totalIonFlux;
volVectorField newIonFlux = 0*ionFlux;
ionDensity = 0*ionDensity;
forAll (speciesList, ii){
	// ============================================== DRIFT DIFFUSION MODEL ==================================================
	if (transportModel[ii] == "driftDiffusion"){
		int isElectron = 0;
		if (ii==0){int isElectron = 0;}
		fvScalarMatrix advEqn
			(   
				fvm::ddt(nDensity[ii])
					+ fvm::div((fvc::interpolate((sign(chargeNumber[ii])*mobility[ii]*EField)) & mesh.Sf()),nDensity[ii],"div(vel,nDensity)") 
					- fvc::laplacian(isElectron*nDensity[0]*mobility[0]*kBoltz/eCharge,eTemp)
					- fvm::laplacian(diffusion[ii], nDensity[ii],"laplacian(difC,nDensity)") + fvm::SuSp(-pSource[ii]/nDensity[ii],nDensity[ii])
			);
		advEqn.solve();
		(nDensity[ii]).max(cosmicRadiation);
		Flux_[ii]=(sign(chargeNumber[ii])*mobility[ii]*EField)*nDensity[ii] - diffusion[ii]*fvc::grad(nDensity[ii])-isElectron*nDensity[0]*mobility[0]*kBoltz/eCharge*fvc::grad(eTemp);
		if (ii==0)
			totalElectronFlux=totalElectronFlux+chargeNumber[ii]*eCharge*Flux_[ii];
		else
			totalIonFlux=totalIonFlux+chargeNumber[ii]*eCharge*Flux_[ii];
			newIonFlux = newIonFlux + Flux_[ii];
		velocity[ii] = Flux_[ii]/(nDensity[ii]+dimensionCorrection*1e-50); //just approx for plot.
	}	
	// ================================================== MOMENTUM MODEL ==================================================
	else if (transportModel[ii] == "momentum"){
		//phie= fvc::interpolate(velocity[ii]) & mesh.Sf();
		fvScalarMatrix advEqn
			(   
				fvm::ddt(nDensity[ii])+
					fvm::div((fvc::interpolate(velocity[ii]) & mesh.Sf()),nDensity[ii],"div(vel,nDensity)") 
					+ fvm::SuSp(-pSource[ii]/nDensity[ii],nDensity[ii])
			);
		advEqn.solve();
		(nDensity[ii]).max(cosmicRadiation);
		Flux_[ii]=velocity[ii]*nDensity[ii];
		if (ii==0)
			totalElectronFlux=totalElectronFlux+chargeNumber[ii]*eCharge*Flux_[ii];
		else
			totalIonFlux=totalIonFlux+chargeNumber[ii]*eCharge*Flux_[ii];
			newIonFlux = newIonFlux + Flux_[ii];
	}
	ionFlux = 1*newIonFlux;
	currentDensity = totalElectronFlux + totalIonFlux ;
	//====================================================================================================
	if (ii > 0)//not electron
		ionDensity = ionDensity + chargeNumber[ii]*nDensity[ii];
	//sigma_impedance = eCharge*eCharge*nDensity[0]*collFreqxN[0]/(9.109e-31*dimensionCorrectionMass*(probingFrequency+collFreqxN[0]/nDensity[0]));
	//sigma_impedance.correctBoundaryConditions();
	//epsilon_impedance = 4*(22/7)*eCharge*nDensity[0]/(9.109e-31*dimensionCorrectionMass*(probingFrequency+collFreqxN[0]/nDensity[0]*collFreqxN[0]/nDensity[0]));
	//epsilon_impedance.correctBoundaryConditions();
}
