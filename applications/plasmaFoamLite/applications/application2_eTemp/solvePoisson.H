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