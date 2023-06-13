/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "eTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eTemp, 0);
    defineRunTimeSelectionTable(eTemp, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eTemp::eTemp
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E
)
:
    thermo_(thermo),
    mspm_(mspm),
    E_(E),
    mesh_(thermo.T().mesh()),
	runTime_(const_cast<Time&>(mesh_.time())),
    restartcapable(runTime_.controlDict().lookup("restartCapable"))
{}


// ************************************************************************* //
