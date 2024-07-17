//*******************************************************************************************************
//================================== ENERGY CONTINUITY ==================================================
//=======================================================================================================
//collisionFreqxN.correctBoundaryConditions();
kappa_e = 3.2*kBoltz*kBoltz*eTemp*nDensity[0]*mobility[0]/eCharge;///(collFreqxN[0]/nDensity[0]*speciesMass[0]*dimensionCorrectionMass); 
kappa_e.correctBoundaryConditions();
jouleHeating= -(EField & Flux_[0])*eCharge;
totEnergyGen=-eSource+jouleHeating;
if (transportModel[0]=="momentum"){//suggested by professor (Derived from Bittencourt book)
	keterm1 = 0.5*nDensity[0]*(velocity[0] & velocity[0])*speciesMass[0]*dimensionCorrectionMass;
	keterm2 = 0.5*nDensity[0]*velocity[0]*(velocity[0]&velocity[0])*speciesMass[0]*dimensionCorrectionMass;
	fvScalarMatrix advEEqn
		(   
			fvm::ddt((3/2)*kBoltz*nDensity[0],eTemp) +
				fvm::div((fvc::interpolate((5/2)*kBoltz*(Flux_[0])) & mesh.Sf()),eTemp,"div(vel,eDensity)") 
				-fvm::laplacian(kappa_e, eTemp,"laplacian(difC,eDensity)") 
				+ fvm::SuSp(-totEnergyGen/eTemp, eTemp) + fvc::ddt(keterm1) + fvc::div(keterm2)
		);
	advEEqn.solve();
	eTemp.max(1e-50);
	//no need to calculate eDensity as it wont be used anywhere in the code.
}
else if (transportModel[0]=="driftDiffusion"){
	fvScalarMatrix advEEqn //energy balance equation identical to Wilson et. al -- 
	(   
		fvm::ddt(eDensity) +
			fvm::div((fvc::interpolate((5/3)*(sign(chargeNumber[0])*mobility[0]*EField)) & mesh.Sf()),eDensity,"div(vel,eDensity)") 
			-fvm::laplacian((5/3)*diffusion[0], eDensity,"laplacian(difC,eDensity)") 
			//-fvc::laplacian(kappa_e, eTemp,"laplacian(difC,eDensity)") 
			+ fvm::SuSp(-totEnergyGen/eDensity, eDensity) 
	);
	advEEqn.solve();
	eDensity.max(1e-50);
	eTemp=eDensity/(3/2)/kBoltz/(nDensity[0]+1e-50*dimensionCorrection);//eTemp used in interpolation
	eTemp.correctBoundaryConditions();
}

//================================= POWER CALCULATION==================================================
//=====================================================================================================
gradeD = fvc::grad(eDensity);
gradne = fvc::grad(nDensity[0]);
powerDissipatedElectron = ((totalElectronFlux) & EField);
powerDissipatedIon = ((totalIonFlux) & EField);
powerDissipated = meshV*((totalElectronFlux + totalIonFlux) & EField);
powerDissipated.correctBoundaryConditions();
power = gSum(powerDissipated)*unitMatrixx*dimensionCorrectionPower;
power.correctBoundaryConditions();