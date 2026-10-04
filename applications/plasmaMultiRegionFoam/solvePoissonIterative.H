{
    emc->correct(voltageDict);

    // Plasma equation of this time step
    volScalarField rhs
    (
        poissonEquationSolver == "semiImplicit"
      ? -eCharge*mspm().potentialExpSource()
      : -eCharge*mspm().netCharge()
    );

    epsEff.internalField() = eps.internalField();

    if (poissonEquationSolver == "semiImplicit")
    {
        epsEff.internalField() +=
            runTime.deltaTValue()*eCharge.value()
           *mspm().potentialImpSource()().internalField();
    }

    epsEff.correctBoundaryConditions();

    // Dirichlet-Neumann iteration: the plasma is solved with the interface
    // potential fixed, each dielectric with the gradient that satisfies
    // Gauss's law, and the interface potential is relaxed towards the one
    // the dielectric obtains
    List<scalarField> oldResidual(couplingPatches.size());

    label iter = 0;
    scalar change = GREAT;

    do
    {
        solve(fvm::laplacian(epsEff, Phi, "laplacian(eps,Phi)") == rhs);

        scalar phiScale = gMax(mag(Phi.internalField()));

        forAll(dielectricRegions, i)
        {
            solve(fvm::laplacian(epsEffD[i], PhiD[i], "laplacian(eps,Phi)"));

            phiScale = max(phiScale, gMax(mag(PhiD[i].internalField())));
        }

        // Residual of the interface potential and Aitken's factor
        List<scalarField> residual(couplingPatches.size());

        scalar num = 0;
        scalar den = 0;

        change = 0;

        forAll(couplingPatches, k)
        {
            const iterativeCoupledPotentialFvPatchScalarField& pPhi =
                refCast<const iterativeCoupledPotentialFvPatchScalarField>
                (
                    Phi.boundaryField()[couplingPatches[k]]
                );

            residual[k] = pPhi.shadowInterfaceValue() - pPhi.refValue();

            change = max(change, gMax(mag(residual[k])));

            if (iter > 0)
            {
                const scalarField dr(residual[k] - oldResidual[k]);

                num += gSum(oldResidual[k]*dr);
                den += gSum(dr*dr);
            }
        }

        if (iter > 0 && den > VSMALL)
        {
            couplingRelaxation *= -num/den;

            couplingRelaxation =
                sign(couplingRelaxation)
               *min(max(mag(couplingRelaxation), 1e-4), 2.0);
        }

        forAll(couplingPatches, k)
        {
            iterativeCoupledPotentialFvPatchScalarField& pPhi =
                refCast<iterativeCoupledPotentialFvPatchScalarField>
                (
                    Phi.boundaryField()[couplingPatches[k]]
                );

            pPhi.refValue() += couplingRelaxation*residual[k];

            oldResidual[k] = residual[k];
        }

        change /= max(phiScale, SMALL);
    }
    while (change > couplingTolerance && ++iter < couplingMaxIterations);

    Info<< "Plasma-dielectric iterations " << min(iter + 1, couplingMaxIterations)
        << ", interface potential change " << change
        << ", relaxation " << couplingRelaxation << endl;

    if (change > couplingTolerance)
    {
        Info<< "Plasma-dielectric iteration not converged to "
            << couplingTolerance << endl;
    }

    Phi.correctBoundaryConditions();

    E = -fvc::grad(Phi);

    if (runTime.outputTime())
    {
        E.write();
    }

    forAll(dielectricRegions, i)
    {
        PhiD[i].correctBoundaryConditions();

        ED[i] = -fvc::grad(PhiD[i]);

        JtotD[i] =
            epsD[i]*(ED[i] - ED[i].oldTime())/runTime.deltaT().value();
    }
}
