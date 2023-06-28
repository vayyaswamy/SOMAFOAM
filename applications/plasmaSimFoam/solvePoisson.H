{
	emc->correct(voltageDict);

	if (poissonEquationSolver == "semiImplicit")
	{
		volScalarField epsSource = epsilon0 + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();

		fvScalarMatrix phiPlasmaEqn
		(
			fvm::laplacian(epsSource, Phi, "laplacian(eps0,Phi)") == -eCharge*(mspm().potentialExpSource())
		);

		phiPlasmaEqn.solve();

		E = -fvc::reconstruct(phiPlasmaEqn.flux())/epsSource;
	}
	else if (poissonEquationSolver == "explicit")
	{
		fvScalarMatrix phiPlasmaEqn
		(
			fvm::laplacian(epsilon0, Phi, "laplacian(eps0,Phi)") == -eCharge*(mspm().netCharge())
		);

		phiPlasmaEqn.solve();

		E = -fvc::reconstruct(phiPlasmaEqn.flux())/epsilon0;
	}
	else if (poissonEquationSolver == "linearized")
	{
		fvScalarMatrix phiPlasmaEqn
		(
			fvm::laplacian(epsilon0, Phi, "laplacian(eps0,Phi)") == -eCharge*(mspm().netCharge() + runTime.deltaTValue()*fvc::div(F.oldTime()))
		);

		phiPlasmaEqn.solve();

		E = -fvc::reconstruct(phiPlasmaEqn.flux())/epsilon0; 
	}
}
