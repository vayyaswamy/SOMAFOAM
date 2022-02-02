/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Namespace
    plasmaConstants

Description
    plasma constants.

\*---------------------------------------------------------------------------*/

#ifndef plasmaConstants_H
#define plasmaConstants_H

#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace plasmaConstants
{
    extern const dimensionedScalar A;

    extern const dimensionedScalar eCharge;

    extern const dimensionedScalar KBE;

    extern const dimensionedScalar boltzC;

    extern const dimensionedScalar boltzCsqr;

    extern const dimensionedScalar rA;

    extern const dimensionedScalar eChargeA;

    extern const dimensionedScalar epsilon0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
