/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "none.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace emcModels
{
    defineTypeNameAndDebug(none, 0);
    addToRunTimeSelectionTable(emcModel, none, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emcModels::none::none
(
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime
)
:
    emcModel(electroMagnetics, mspm, E, runTime)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::emcModels::none::~none()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::emcModels::none::correct(dictionary& voltageDict)
{}

bool Foam::emcModels::none::read(const dictionary& electroMagnetics)
{
    emcModel::read(electroMagnetics);

    return true;
}
// ************************************************************************* //
