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

while (pimple.correct())
{
        #include "solvePoissonD.H"
}
#include "detachPatches.H"

#include "plasmaEqn.H"

#include "surfaceCharge.H"
