/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Description
	plasma model constants
\*---------------------------------------------------------------------------*/

#include "plasmaConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Avogadro constant (number of molecules in a mole)(gram corrected)
const Foam::dimensionedScalar Foam::plasmaConstants::A
(
    "A",
    dimensionSet(0, 0, 0, 0, -1, 0, 0),
    6.02214076e26
);

// Elementery charge (charge of electron)
const Foam::dimensionedScalar Foam::plasmaConstants::eCharge
(
    "eCharge",
    dimensionSet(0, 0, 1, 0, 0, 1, 0),
    1.602176487e-19
);

// KB/eCharge
const Foam::dimensionedScalar Foam::plasmaConstants::KBE
(
    "KBE",
    dimensionSet(1, 2, -3, -1, 0, -1, 0),
    8.61733034152e-05
);

const Foam::dimensionedScalar Foam::plasmaConstants::boltzC
(
    "boltzC",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    1.38064852e-23
);

const Foam::dimensionedScalar Foam::plasmaConstants::boltzCsqr
(
    "boltzCsqr",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    1.90619033577819e-46
);

//Avogadro constant reciprocal(gram corrected)
const Foam::dimensionedScalar Foam::plasmaConstants::rA
(
    "rA",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    1.6605392767355e-27
);

//eCharge*A
const Foam::dimensionedScalar Foam::plasmaConstants::eChargeA
(
    "eChargeA",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    96485323.27076
);

//permittivity of free space
const Foam::dimensionedScalar Foam::plasmaConstants::epsilon0
(
    "epsilon0",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    8.85418782E-12
);

// ************************************************************************* //
