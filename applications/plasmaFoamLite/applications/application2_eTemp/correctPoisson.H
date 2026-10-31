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