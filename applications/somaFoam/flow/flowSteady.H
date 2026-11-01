// Steady gas flow by the SIMPLE iterations of rhoSimpleFoam, with the plasma
// force and heating of plasmaForcePtr and plasmaHeatPtr, until U, p (and T)
// change by less than steadyTolerance in an iteration, or for
// steadyMaxIterations iterations
{
    volVectorField& U = UgasPtr();
    surfaceScalarField& phi = phiPtr();
    volScalarField& p = thermo.p();
    volScalarField& hs = thermo.hs();
    const volScalarField& psi = thermo.psi();
    const volVectorField& plasmaForce = plasmaForcePtr();
    const volScalarField& plasmaHeat = plasmaHeatPtr();

    scalar flowChange = GREAT;
    label iter = 0;

    while (iter < steadyMaxIterations && flowChange >= steadyTolerance)
    {
        iter++;

        const volVectorField U0("U0", U);
        const volScalarField p0("p0", p);
        const volScalarField T0("T0", thermo.T());

        // The relaxation of p and rho uses their previous iteration
        p.storePrevIter();
        rho.storePrevIter();

        #include "steadyUEqn.H"

        if (flowEnergy)
        {
            #include "steadyHsEqn.H"
        }

        #include "steadyPEqn.H"

        turbulence->correct();

        flowChange = max
        (
            gMax(mag(U.internalField() - U0.internalField()))
           /max(gMax(mag(U.internalField())), SMALL),
            gMax(mag(p.internalField() - p0.internalField()))
           /max(gMax(mag(p.internalField())), SMALL)
        );

        if (flowEnergy)
        {
            flowChange = max
            (
                flowChange,
                gMax(mag(thermo.T().internalField() - T0.internalField()))
               /gMax(thermo.T().internalField())
            );
        }
    }

    // Volume flux that carries the species
    phiGasPtr() = phi/linearInterpolate(rho);

    Info<< "Steady flow: " << iter << " SIMPLE iterations, last change "
        << flowChange << (flowChange < steadyTolerance ? "" : " (NOT converged)")
        << ", max |U| = " << gMax(mag(U.internalField()))
        << " m/s, T = " << gMin(thermo.T().internalField()) << " to "
        << gMax(thermo.T().internalField()) << " K" << endl;
}
