{
    emc->correct(voltageDict);

    if (poissonEquationSolver == "semiImplicit")
    {
        volScalarField epsSource = eps + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();


        fvScalarMatrix phiPlasmaEqn
        (
            fvm::laplacian(epsSource, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().potentialExpSource())
        );


        phiPlasmaEqn.solve();


        E = -fvc::grad(Phi);

    }
    else if (poissonEquationSolver == "explicit")
    {
        fvScalarMatrix phiPlasmaEqn
        (
            fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().netCharge())
        );

        phiPlasmaEqn.relax();

        phiPlasmaEqn.solve();

        relaxFinal(Phi);


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
