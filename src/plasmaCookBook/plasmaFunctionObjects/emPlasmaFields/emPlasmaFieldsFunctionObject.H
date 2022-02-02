/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Typedef
    Foam::emPlasmaFieldsFunctionObject

Description
    FunctionObject wrapper around emPlasmaFields to allow it to be created via
    the functions list within controlDict.

SourceFiles
    emPlasmaFieldsFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef emPlasmaFieldsFunctionObject_H
#define emPlasmaFieldsFunctionObject_H

#include "emPlasmaFields.H"
#include "OutputFilterFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef OutputFilterFunctionObject<emPlasmaFields>
        emPlasmaFieldsFunctionObject;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
