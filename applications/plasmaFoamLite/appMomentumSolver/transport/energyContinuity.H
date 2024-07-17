//*******************************************************************************************************
//================================== ENERGY CONTINUITY ==================================================
//=======================================================================================================
//collisionFreqxN.correctBoundaryConditions();
gradTe = fvc::grad(eTemp);
kappa_e = 3.2*kBoltz*kBoltz*eTemp*nDensity[0]*mobility[0]/eCharge;///(collFreqxN[0]/nDensity[0]*speciesMass[0]*dimensionCorrectionMass); 
kappa_e.correctBoundaryConditions();
jouleHeating= -(EField & Flux_[0])*eCharge;
totEnergyGen=-eSource+jouleHeating;
int isMomentum = 0;
if (transportModel[0]=="momentum")
	isMomentum = 1 ;

//suggested by professor (Derived from Bittencourt book)
keterm1 = 0.5*nDensity[0]*(velocity[0] & velocity[0])*speciesMass[0]*dimensionCorrectionMass;
keterm2 = 0.5*nDensity[0]*velocity[0]*(velocity[0]&velocity[0])*speciesMass[0]*dimensionCorrectionMass;
fvScalarMatrix advEEqn
	(   
		fvm::ddt((3/2)*kBoltz*nDensity[0],eTemp) +
			fvm::div((fvc::interpolate((5/2)*kBoltz*(Flux_[0])) & mesh.Sf()),eTemp,"div(vel,eDensity)") 
			-fvm::laplacian(kappa_e, eTemp,"laplacian(difC,eDensity)") 
			+ fvm::SuSp(-totEnergyGen/(eTemp+1e-50*dimensionCorrectionTemp), eTemp) + isMomentum*(fvc::ddt(keterm1) + fvc::div(keterm2))
	);
advEEqn.solve();
eTemp.max(300);

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
