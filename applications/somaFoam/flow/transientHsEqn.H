// Energy equation of rhoPimpleFoam (foam-extend 4.0) for the sensible
// enthalpy hs of the plasma thermo, with the heating of the plasma and a
// uniform heating (flowProperties heatSource)
{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + fvm::div(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        DpDtPtr()
      + plasmaHeat
      + flowHeatSource
    );

    hsEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor(hs.select(flowPimple().finalIter()))
    );

    hsEqn.solve(mesh.solutionDict().solver((hs.select(flowPimple().finalIter()))));

    thermo.correct();
}
