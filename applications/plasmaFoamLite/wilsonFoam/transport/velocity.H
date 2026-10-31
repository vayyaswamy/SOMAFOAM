//*****************************************************************************************************
//====================================== MOMENTUM SOLVER ==============================================
//*****************************************************************************************************
//ionFlux=0*ionFlux;
forAll (speciesList, ii){
	if (collisionType[ii] == "muBased" && chargeNumber[ii]!=0) //muBased done here. rrBased is default
		collFreqxN[ii] = nDensity[ii]*eCharge/(mobility[ii])/speciesMass[ii]/dimensionCorrectionMass;
	if (transportModel[ii] == "momentum"){
		if (ii==0){//electron
			pTemperature= 1*eTemp;}
		else //ion and other species
			{pTemperature= gasTemperature*dimensionCorrectionTemp;}
		Info << "Solving momentum equation" << endl;//--------------------
			fvVectorMatrix momEqn
			(
				fvm::ddt(nDensity[ii],velocity[ii])  
				+ fvm::div(fvc::interpolate(nDensity[ii])*(fvc::interpolate(velocity[ii]) & mesh.Sf()),velocity[ii],"div(velocity*nDensity,velocity)")  
				+ fvm::Sp(collFreqxN[ii],velocity[ii])
				== 
				chargeNumber[ii]/abs(chargeNumber[ii])*eCharge*nDensity[ii]*EField/speciesMass[ii]/dimensionCorrectionMass
				- fvc::grad(kBoltz*pTemperature*nDensity[ii])/speciesMass[ii]/dimensionCorrectionMass
			);
		momEqn.solve();
	}
}//======================================================================================================
eCollisionFrequency = collFreqxN[0]/nDensity[0];