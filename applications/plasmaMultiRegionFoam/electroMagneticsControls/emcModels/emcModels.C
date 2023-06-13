/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "emcModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(emcModel, 0);
    defineRunTimeSelectionTable(emcModel, dictionary);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emcModel::emcModel
(
    const dictionary& electroMagnetics,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E,
	const Time& runTime
)
:
    regIOobject
    (
        IOobject
        (
            "emcModel",
            E.time().constant(),
            E.db()
        )
    ),
    emcModelCoeffs_
    (
        electroMagnetics.subDict
        (
            word(electroMagnetics.lookup("emcModel")) + "Coeffs"
        )
    ),
    mspm_(mspm),
    E_(E),
	time_(runTime)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::emcModel::~emcModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::emcModel::read(const dictionary& electroMagnetics)
{
    emcModelCoeffs_ = electroMagnetics.subDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
