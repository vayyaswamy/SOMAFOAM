//=================== semi implicit poisson solver ================================
if (emSolver=="semiImplicit"){
	netCharge=0*netCharge;
	impSource= 0*impSource;
	forAll(speciesList, ii){
					netCharge = netCharge + eCharge*chargeNumber[ii]*nDensity[ii] + 
							eCharge*chargeNumber[ii]*runTime.deltaT()*fvc::laplacian(diffusion[ii], nDensity[ii],"laplacian(difC,nDensity)");
					impSource = impSource +  eCharge*chargeNumber[ii]*nDensity[ii]*mobility[ii]*runTime.deltaT()/epsilon;
	}
	impSource=impSource+epsilon/epsilon;

	fvScalarMatrix poissonPlasmaEquation
							(
									fvm::laplacian(impSource,EPotential) + netCharge/epsilon
							);
	poissonPlasmaEquation.relax();
	poissonPlasmaEquation.solve();
	EField=-fvc::grad(EPotential);
}//--------------------------------------------------------------------------------

// =================== Explicit poisson solver ====================================
else{
	netCharge=0*netCharge;
	forAll(speciesList, ii){
			netCharge = netCharge + eCharge*chargeNumber[ii]*nDensity[ii];
	}
	fvScalarMatrix poissonPlasmaEquation
					(
						fvm::laplacian(EPotential) == -netCharge/epsilon
					);
	poissonPlasmaEquation.relax();
	poissonPlasmaEquation.solve();
	EField=-fvc::grad(EPotential);
	Info << "solvedPoisson" << endl;
}//--------------------------------------------------------------------------------