// Pressure equation of rhoPimpleFoam (foam-extend 4.0); the PIMPLE controls
// are those of flowPimple (flowPIMPLE in fvSolution)

rho = thermo.rho();

volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();

if (flowPimple().nCorrPISO() <= 1)
{
    UEqn.clear();
}

if (flowPimple().transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    while (flowPimple().correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve
        (
            mesh.solutionDict().solver(p.select(flowPimple().finalInnerIter()))
        );

        if (flowPimple().finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
        //+ fvc::ddtPhiCorr(rUA, rho, U, phi)
        );

    while (flowPimple().correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve
        (
            mesh.solutionDict().solver(p.select(flowPimple().finalInnerIter()))
        );

        if (flowPimple().finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

//if (oCorr != nOuterCorr-1)
{
    // Explicitly relax pressure for momentum corrector
    p.relax();

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value()
        << " " << min(rho).value() << endl;
}

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

DpDtPtr() = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

bound(p, pMin);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
/*
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}
*/
