/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Typedef
    Foam::IOemPlasmaFields

Description

\*---------------------------------------------------------------------------*/

#ifndef IOemPlasmaFields_H
#define IOemPlasmaFields_H

#include "emPlasmaFields.H"
#include "IOOutputFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef IOOutputFilter<emPlasmaFields> IOemPlasmaFields;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
