/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "makeMultiSpeciesPlasmaModel.H"
#include "driftDiffusion.H"
#include "driftDiffusionOS.H"
#include "mixed.H"
//#include "mixedN.H"
#include "momentum.H"
#include "zeroD.H"
#include "thermoPhysicsTypes.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makeMultiSpeciesPlasmaModel(driftDiffusion, constGasThermoPhysics);

makeMultiSpeciesPlasmaModel(momentum, constGasThermoPhysics);

makeMultiSpeciesPlasmaModel(mixed, constGasThermoPhysics);

//makeMultiSpeciesPlasmaModel(mixedN, constGasThermoPhysics);

makeMultiSpeciesPlasmaModel(zeroD, constGasThermoPhysics);

makeMultiSpeciesPlasmaModel(driftDiffusionOS, constGasThermoPhysics);
}

// ************************************************************************* //
