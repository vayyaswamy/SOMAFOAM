// One time step of rhoPimpleFoam (the plasma step), with the plasma force and
// heating of plasmaForcePtr and plasmaHeatPtr
{
    volVectorField& U = UgasPtr();
    surfaceScalarField& phi = phiPtr();
    volScalarField& p = thermo.p();
    volScalarField& hs = thermo.hs();
    const volScalarField& psi = thermo.psi();
    const volVectorField& plasmaForce = plasmaForcePtr();
    const volScalarField& plasmaHeat = plasmaHeatPtr();

    #include "rhoEqn.H"

    // --- Pressure-velocity PIMPLE corrector loop
    while (flowPimple().loop())
    {
        #include "transientUEqn.H"

        if (flowEnergy)
        {
            #include "transientHsEqn.H"
        }

        // --- PISO loop
        while (flowPimple().correct())
        {
            #include "transientPEqn.H"
        }

        turbulence->correct();
    }

    // Volume flux that carries the species
    phiGasPtr() = phi/linearInterpolate(rho);

    Info<< "Flow: max |U| = " << gMax(mag(U.internalField())) << " m/s, p = "
        << gMin(p.internalField()) << " to " << gMax(p.internalField())
        << " Pa, T = " << gMin(thermo.T().internalField()) << " to "
        << gMax(thermo.T().internalField()) << " K" << endl;
}
