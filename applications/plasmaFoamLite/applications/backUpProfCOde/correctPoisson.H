// ------------------ NET CHARGE CALCULATION -------------------------------------
netCharge=0*netCharge;
forAll(speciesList, ii){
		netCharge= netCharge + eCharge*chargeNumber[ii]*nDensity[ii];
}
fvScalarMatrix poissonPlasmaEquation
				(
					fvm::laplacian(EPotential) + (netCharge) /epsilon //explicit term
				);
poissonPlasmaEquation.relax();
poissonPlasmaEquation.solve();
EField=-fvc::grad(EPotential);


/*
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
*/
