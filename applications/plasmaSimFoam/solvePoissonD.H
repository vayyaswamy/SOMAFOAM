{
	emc->correct(voltageDict);

	if (poissonEquationSolver == "semiImplicit")
	{
		volScalarField epsSource = epsilon0 + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();

		coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);

		fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
		(
			fvm::laplacian(epsSource, Phi, "laplacian(eps0,Phi)") == -eCharge*(mspm().potentialExpSource()) - surfC
		);

		phiEqns.set(0, phiPlasmaEqn);

		forAll(dielectricRegions, i)
		{
			phiEqns.set(i+1, new fvScalarMatrix
			(
				fvm::laplacian(eps[i], PhiD[i], "laplacian(eps,Phi)") 
			));
		}

		phiEqns.solve();

		E = -fvc::reconstruct(phiPlasmaEqn->flux())/epsSource;

		forAll(dielectricRegions, i)
		{
			ED[i] = -fvc::grad(PhiD[i]);
		}   
	}
	else if (poissonEquationSolver == "explicit")
	{
		coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);

		fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
		(
			fvm::laplacian(epsilon0, Phi, "laplacian(eps0,Phi)") == -(eCharge*mspm().netCharge()) - surfC
		);

		phiEqns.set(0, phiPlasmaEqn);

		forAll(dielectricRegions, i)
		{
			phiEqns.set(i+1, new fvScalarMatrix
			(
				fvm::laplacian(eps[i], PhiD[i], "laplacian(eps,Phi)") 
			));
		}

		phiEqns.solve();

		E = -fvc::reconstruct(phiPlasmaEqn->flux())/epsilon0;

		forAll(dielectricRegions, i)
		{
			ED[i] = -fvc::grad(PhiD[i]);
		}   
	}
	else if (poissonEquationSolver == "linearized")
	{
		coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);

		fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
		(
			fvm::laplacian(epsilon0, Phi, "laplacian(eps0,Phi)") == -eCharge*(mspm().netCharge() + runTime.deltaTValue()*fvc::div(F.oldTime())) - surfC
		);

		phiEqns.set(0, phiPlasmaEqn);

		forAll(dielectricRegions, i)
		{
			phiEqns.set(i+1, new fvScalarMatrix
			(
				fvm::laplacian(eps[i], PhiD[i], "laplacian(eps,Phi)") 
			));
		}

		phiEqns.solve();

		E = -fvc::reconstruct(phiPlasmaEqn->flux())/epsilon0;

		forAll(dielectricRegions, i)
		{
			ED[i] = -fvc::grad(PhiD[i]);
		}   
	}
}
