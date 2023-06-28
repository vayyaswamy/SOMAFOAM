/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(plasmaEnergyModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaEnergyModel::plasmaEnergyModel
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E
)
:
    IOdictionary
    (
        IOobject
        (
            "plasmaProperties",
			thermo.T().mesh().time().constant(),
			thermo.T().mesh(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
        )
    ),
    eTempPtr_(eTemp::New(thermo, mspm, E, subDict("energyModel"))),
    iTempPtr_(iTemp::New(thermo, mspm, E, subDict("energyModel"))),
    gTempPtr_(gTemp::New(thermo, mspm, E, subDict("energyModel")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar plasmaEnergyModel::ecorrect(psiChemistryModel& chemistry, const volVectorField& E)
{
    scalar output = eTempPtr_->correct(chemistry, E);
    return output;
}

void plasmaEnergyModel::icorrect(psiChemistryModel& chemistry)
{
    iTempPtr_->correct(chemistry);
}

void plasmaEnergyModel::gcorrect(psiChemistryModel& chemistry, const volVectorField& E)
{
    gTempPtr_->correct(chemistry, E);
}

bool plasmaEnergyModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
