{
	emc->correct(voltageDict);

	if (poissonEquationSolver == "semiImplicit")
	{
		volScalarField epsSource = eps + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();

		//Info << "Inside Poisson" << endl;

		fvScalarMatrix phiPlasmaEqn
		(
			fvm::laplacian(epsSource, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().potentialExpSource())
		);

		//Info << "About to solve " << endl;

		phiPlasmaEqn.solve();

		//Info << "Solve done " << endl;

		//E = -fvc::reconstruct(phiPlasmaEqn.flux())/epsSource;
		E = -fvc::grad(Phi);

		//Info << "E-field done" << endl;
	}
	else if (poissonEquationSolver == "explicit")
	{
		fvScalarMatrix phiPlasmaEqn
		(
			fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().netCharge())
		);

		phiPlasmaEqn.relax();

		phiPlasmaEqn.solve();

		Phi.relax();

		//E = -fvc::reconstruct(phiPlasmaEqn.flux())/eps;

		E = -fvc::grad(Phi);
	}
	else if (poissonEquationSolver == "linearized")
	{
		fvScalarMatrix phiPlasmaEqn
		(
			fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().netCharge())//+ runTime.deltaTValue()*fvc::div(F.oldTime()))
		);

		phiPlasmaEqn.solve();

		E = -fvc::reconstruct(phiPlasmaEqn.flux())/eps; 
	}
}
