{
    emc->correct(voltageDict);


    scalar relativePhiResidual = 1.0;

    int iCorr = 0;


    do
    {

    Phi.storePrevIter();

    if (poissonEquationSolver == "semiImplicit")
    {


        coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);


        eps.correctBoundaryConditions();

        volScalarField epsSource = eps + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();


        fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
        (
            fvm::laplacian(epsSource, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().potentialExpSource())

        );


        phiEqns.set(0, phiPlasmaEqn);


        forAll(dielectricRegions, i)
        {
            epsD[i].correctBoundaryConditions();

        }


        forAll(dielectricRegions, i)
        {

            phiEqns.set(i+1, new fvScalarMatrix
            (
                fvm::laplacian(epsD[i], PhiD[i], "laplacian(eps,Phi)")
            ));
        }


        phiEqns.solve();


    }
    else if (poissonEquationSolver == "explicit")
    {
        coupledFvScalarMatrix phiEqns(dielectricRegions.size()+1);


        eps.correctBoundaryConditions();

        forAll(dielectricRegions, i)
        {
            epsD[i].correctBoundaryConditions();

        }


        fvScalarMatrix* phiPlasmaEqn = new fvScalarMatrix
        (
            fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -(eCharge*mspm().netCharge())
        );

        phiEqns.set(0, phiPlasmaEqn);


        forAll(dielectricRegions, i)
        {

            phiEqns.set(i+1, new fvScalarMatrix
            (
                fvm::laplacian(epsD[i], PhiD[i], "laplacian(eps,Phi)")
            ));
        }


        phiEqns.solve();


        E = -fvc::grad(Phi);

        forAll(dielectricRegions, i)
        {
            ED[i] = -fvc::grad(PhiD[i]);
        }

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

    Info << "Calculating E " << endl;

    Phi.correctBoundaryConditions();

    E = -fvc::grad(Phi);

    E.write();

    Info << "Calculating E done " << endl;


    forAll(dielectricRegions, i)
    {
        ED[i] = -fvc::grad(PhiD[i]);
    }

    forAll(dielectricRegions, i)
    {

        JtotD[i] = epsD[i]*(ED[i] - ED[i].oldTime())/runTime.deltaT().value() ;

    }

}
