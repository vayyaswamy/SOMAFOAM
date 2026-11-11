{
    emc->correct(voltageDict);

    if (poissonEquationSolver == "none")
    {
        // The plasma model sets the field (zeroD with a zeroDDrive)
        mspm().setDrivenField(E);
    }
    else if (poissonEquationSolver == "semiImplicit")
    {
        volScalarField epsSource = eps + runTime.deltaTValue()*eCharge*mspm().potentialImpSource();


        for (label nonOrth = 0; nonOrth <= mspm().nNonOrthCorr(); nonOrth++)
        {
            fvScalarMatrix phiPlasmaEqn
            (
                fvm::laplacian(epsSource, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().potentialExpSource())
            );

            phiPlasmaEqn.solve();
        }


        E = -fvc::grad(Phi);

    }
    else if (poissonEquationSolver == "explicit")
    {
        Phi.storePrevIter();

        for (label nonOrth = 0; nonOrth <= mspm().nNonOrthCorr(); nonOrth++)
        {
            fvScalarMatrix phiPlasmaEqn
            (
                fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().netCharge())
            );

            phiPlasmaEqn.relax();

            phiPlasmaEqn.solve();
        }

        relaxFinal(Phi);


        E = -fvc::grad(Phi);
    }
    else if (poissonEquationSolver == "linearized")
    {
        for (label nonOrth = 0; nonOrth <= mspm().nNonOrthCorr(); nonOrth++)
        {
            fvScalarMatrix phiPlasmaEqn
            (
                fvm::laplacian(eps, Phi, "laplacian(eps,Phi)") == -eCharge*(mspm().netCharge())//+ runTime.deltaTValue()*fvc::div(F.oldTime()))
            );

            phiPlasmaEqn.solve();

            if (nonOrth == mspm().nNonOrthCorr())
            {
                E = -fvc::reconstruct(phiPlasmaEqn.flux())/eps;
            }
        }
    }
}
