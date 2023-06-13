/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "gasCorrectedTemp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gasCorrectedTemp, 0);
    addToRunTimeSelectionTable(iTemp, gasCorrectedTemp, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasCorrectedTemp::gasCorrectedTemp
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E,
    const dictionary& dict
)
:
    iTemp(thermo, mspm, E),
    ionSpecie(dict.lookup("ionSpecie")),
    iIndex_(mspm.species()[ionSpecie]),
    backgroundGas(dict.lookup("backgroundGas")),
    bIndex_(mspm.species()[backgroundGas])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gasCorrectedTemp::correct
(
    psiChemistryModel& chemistry
)
{
	volScalarField& TiC = thermo().Tion();

	const volScalarField& TgC = thermo().T();

    TiC = TgC + (pow((mspm().mu(iIndex_)*Foam::mag(E())),2)*(plasmaConstants::rA*mspm().W(bIndex_)/plasmaConstants::boltzC)*(mspm().W(bIndex_)+mspm().W(iIndex_)/(3*mspm().W(bIndex_)+5*mspm().W(iIndex_))));
}


// ************************************************************************* //
