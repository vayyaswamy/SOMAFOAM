// Momentum equation of rhoPimpleFoam (foam-extend 4.0), with the force of the
// plasma and a uniform force (flowProperties bodyForce); the PIMPLE controls
// are those of flowPimple (flowPIMPLE in fvSolution)

// Solve the Momentum equation
tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff()
 ==
    plasmaForce
  + flowBodyForce
);

UEqn().relax
(
    mesh.solutionDict().equationRelaxationFactor(U.select(flowPimple().finalIter()))
);

volScalarField rUA = 1.0/UEqn().A();

if (flowPimple().momentumPredictor())
{
    solve
    (
        UEqn() == -fvc::grad(p),
        mesh.solutionDict().solver((U.select(flowPimple().finalIter())))
    );
}
else
{
    U = rUA*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}
