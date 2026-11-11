// One outer iteration of the plasma-dielectric step: Poisson equation with
// the dielectrics attached, then the plasma equations and the surface charge
#include "attachPatches.H"

if (dielectricWeightsOutdated)
{
    forAll(dielectricRegions, i)
    {
        dielectricRegions[i].surfaceInterpolation::movePoints();

        // Distance between the cell centres across the
        // interface, used by the harmonic interpolation of
        // the permittivity (cached per mesh)
        magLongDelta::Delete(dielectricRegions[i]);
    }

    dielectricWeightsOutdated = false;
}

solverTimers::outer++;

double tTimer = solverTimers::now();

while (pimple.correct())
{
        #include "solvePoissonD.H"
}
#include "detachPatches.H"

solverTimers::poisson += solverTimers::now() - tTimer;
tTimer = solverTimers::now();

#include "plasmaEqn.H"

solverTimers::species += solverTimers::now() - tTimer;
tTimer = solverTimers::now();

#include "surfaceCharge.H"

solverTimers::surface += solverTimers::now() - tTimer;
