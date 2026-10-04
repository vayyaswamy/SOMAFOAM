/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Description
    Sources shared with somaFoam (voltage and power control models, region
    properties), compiled here from the somaFoam directory so that the two
    solvers cannot drift apart.

\*---------------------------------------------------------------------------*/

#include "electroMagneticsControls/emcModels/emcModels.C"
#include "electroMagneticsControls/emcModels/newEmcModels.C"
#include "electroMagneticsControls/voltageControl/voltageControl.C"
#include "electroMagneticsControls/powerControl/powerControl.C"
#include "electroMagneticsControls/none/none.C"
#include "regionProperties/regionProperties.C"

// ************************************************************************* //
