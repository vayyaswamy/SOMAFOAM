/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "gTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gTemp, 0);
    defineRunTimeSelectionTable(gTemp, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gTemp::gTemp
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
	runTime_(const_cast<Time&>(mesh_.time()))
{}


// ************************************************************************* //