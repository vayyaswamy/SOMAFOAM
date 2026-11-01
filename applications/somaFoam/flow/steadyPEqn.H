// Pressure equation of rhoSimpleFoam (foam-extend 4.0); the SIMPLE controls
// are those of flowSimple (SIMPLE in fvSolution)

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();
UEqn.clear();

bool closedVolume = false;

if (flowSimple().transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)*(fvc::interpolate(U) & mesh.Sf())
    );

    while (flowSimple().correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (flowSimple().finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    closedVolume = adjustPhi(phi, U, p);

    while (flowSimple().correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (flowSimple().finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }
}


// (continuity errors of the incompressible form not reported)

// Explicitly relax pressure for momentum corrector
p.relax();

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
