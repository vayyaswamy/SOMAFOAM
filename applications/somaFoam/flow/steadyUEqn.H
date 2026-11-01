// Momentum equation of rhoSimpleFoam (foam-extend 4.0), with the force of the
// plasma and a uniform force (flowProperties bodyForce)

    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff()
     ==
        plasmaForce
      + flowBodyForce
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
