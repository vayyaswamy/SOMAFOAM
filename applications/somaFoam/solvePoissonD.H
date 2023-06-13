{
	emc->correct(voltageDict);

	//Info << "solvePoissonD" << endl;

	scalar relativePhiResidual = 1.0;

	int iCorr = 0;



	do
	{

	Phi.storePrevIter();

	if (poissonEquationSolver == "semiImplicit")
	{
		

		coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);

		//Info << "eps plasma " << endl;

		eps.correctBoundaryConditions();

		volScalarField epsSource = eps + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();

		//Info << "eps = " << eps << endl;

		

		//Info << "eps dielectric " << endl;

		
		//Info << "Done with eps" << endl;

		//Info << "eps corrected" << endl;

		//Info << "setting matrix " << endl;

		fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
		(
			fvm::laplacian(epsSource, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().potentialExpSource())

			//fvm::laplacian(epsSource, Phi, "laplacian(eps,Phi)") == 0
		);

		//Info << "intermediate " << endl;

		phiEqns.set(0, phiPlasmaEqn);

		//Info << "Done setting " << endl;

		//Info << "setting dielectric matrix " << endl;

		forAll(dielectricRegions, i)
		{
			epsD[i].correctBoundaryConditions();

		}

		//Info << "epsD = " << epsD << endl;



		forAll(dielectricRegions, i)
		{

			phiEqns.set(i+1, new fvScalarMatrix
			(
				fvm::laplacian(epsD[i], PhiD[i], "laplacian(eps,Phi)") 
			));
		}
		//Info << "Done setting dielectric matrix " << endl;



		phiEqns.solve();

	
	}
	else if (poissonEquationSolver == "explicit")
	{
		coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);

		//Info << "eps = " << eps << endl;

		//Info << "eps correct " << endl;

		eps.correctBoundaryConditions();

		forAll(dielectricRegions, i)
		{
			epsD[i].correctBoundaryConditions();

		}

		//Info << "eps correct done " << endl;

		//Info << "Stage 2 " << endl;

		fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
		(
			fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -(eCharge*mspm().netCharge())
		);

		phiEqns.set(0, phiPlasmaEqn);

		//Info << "Stage 3" << endl;

		//Info << "eps = " << eps << endl;

		//Info << "epsD before correctBoundaryConditions = " << epsD << endl;

		forAll(dielectricRegions, i)
		{

			phiEqns.set(i+1, new fvScalarMatrix
			(
				fvm::laplacian(epsD[i], PhiD[i], "laplacian(eps,Phi)")  
			));
		}

		//Info << "Stage 4" << endl;

		//Info << "epsD = " << epsD << endl;

		//Info << "phiEqns solving " << endl;

		phiEqns.solve();

		//Info << "Done " << endl;

		//Info << "Stage 4a" << endl;

		//Info << "eps = " << eps.boundaryField() << endl;

		//E = -fvc::reconstruct(phiPlasmaEqn->flux())/(eps+SMALL);

		E = -fvc::grad(Phi);

		forAll(dielectricRegions, i)
		{
			ED[i] = -fvc::grad(PhiD[i]);
		}   

		//Info << "Stage 5" << endl;
	}
	else if (poissonEquationSolver == "linearized")
	{
		coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);

		eps.correctBoundaryConditions();

		fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
		(
			fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().netCharge()) //+ runTime.deltaTValue()*fvc::div(F.oldTime()))
		);

		phiEqns.set(0, phiPlasmaEqn);

		forAll(dielectricRegions, i)
		{
			epsD[i].correctBoundaryConditions();

			phiEqns.set(i+1, new fvScalarMatrix
			(
				fvm::laplacian(epsD[i], PhiD[i], "laplacian(eps,Phi)") 
			));
		}

		phiEqns.solve();

		E = -fvc::reconstruct(phiPlasmaEqn->flux())/eps;

		forAll(dielectricRegions, i)
		{
			ED[i] = -fvc::grad(PhiD[i]);
		}   
	}

	#include "checkPhiConvergence.H"

	if (iCorr > 10)
	{

		Info << "iCorr = " << iCorr << "relativePhiResidual = " << relativePhiResidual << endl;	
	}

	

	} 
	while 
	(
	(relativePhiResidual >= 1e-4) 
	&& ++iCorr <= 100
	);

	E = -fvc::grad(Phi);

	//Info << "ddt of E = " << fvc::ddt(E) << endl;

	forAll(dielectricRegions, i)
	{
		ED[i] = -fvc::grad(PhiD[i]);
	}

	forAll(dielectricRegions, i)
	{

		JtotD[i] = epsD[i]*(ED[i] - ED[i].oldTime())/runTime.deltaT().value() ;
		//Info << "ddt of ED in somaFoam = " << JtotD[i] << endl;
			
	}  

}
