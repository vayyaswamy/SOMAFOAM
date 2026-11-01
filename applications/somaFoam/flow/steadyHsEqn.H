// Energy equation of rhoSimpleFoam (foam-extend 4.0) for the sensible
// enthalpy hs of the plasma thermo, with the heating of the plasma and a
// uniform heating (flowProperties heatSource)
{
    fvScalarMatrix hsEqn
    (
        fvm::div(phi, hs)
      - fvm::Sp(fvc::div(phi), hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        fvc::div(phi/fvc::interpolate(rho), rho/psi, "div(U,p)")
      - (rho/psi)*fvc::div(phi/fvc::interpolate(rho))
      + plasmaHeat
      + flowHeatSource
    );

    hsEqn.relax();

    hsEqn.solve();

    thermo.correct();
}
