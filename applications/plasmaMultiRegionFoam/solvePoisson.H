{
	emc->correct(voltageDict);

    int iCorr = 0;

    lduSolverPerformance solverPerf;

    scalar residual = 1.0;

    do
    {
        Phi.storePrevIter();

		if (poissonEquationSolver == "semiImplicit")
		{
			volScalarField epsSource = eps + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();

			fvScalarMatrix phiPlasmaEqn
			(
				fvm::laplacian(epsSource, Phi, "laplacian(eps,phi)")
				==
				- eCharge*(mspm().potentialExpSource())
			);

			phiPlasmaEqn.solve();

			E = - fvc::reconstruct(phiPlasmaEqn.flux())/epsSource;
		}
		else if (poissonEquationSolver == "explicit")
		{
			fvScalarMatrix phiPlasmaEqn
			(
				fvm::laplacian(eps, Phi, "laplacian(eps,phi)")
				==
				- eCharge*(mspm().netCharge())
			);

			phiPlasmaEqn.solve();

			E = -fvc::reconstruct(phiPlasmaEqn.flux())/eps;
		}
		else if (poissonEquationSolver == "linearized")
		{
			fvScalarMatrix phiPlasmaEqn
			(
				fvm::laplacian(eps, Phi, "laplacian(eps,phi)") 
				==
				- eCharge*(mspm().netCharge()
				+ runTime.deltaTValue()*fvc::div(F.oldTime()))
			);

			phiPlasmaEqn.solve();

			E = - fvc::reconstruct(phiPlasmaEqn.flux())/eps;
		}

        // if (iCorr == 0)
        // {
        //     initialResidual = solverPerf.initialResidual();
        // }

		forAll(dielectricRegions, i)
		{
			fvScalarMatrix phiDielectricEqn
			(
				fvm::laplacian(epsD[i], PhiD[i], "laplacian(eps,phi)")
			);

			phiDielectricEqn.solve();

			ED[i] = - fvc::grad(PhiD[i], "grad(phi)");
		}

		residual =
		    gMax
		    (
		        mag
		        (
		            Phi.internalField()
		          - Phi.prevIter().internalField()
		        )
		       /Phi.internalField()
		    );
    }
    while
    (
		mag(residual) > convergenceTolerance && ++iCorr < nCorr
    );

	Info<< "Solving for phi "
	<< ", outer iterations " << iCorr
	<< ", Relative error: " << residual << endl;
}
