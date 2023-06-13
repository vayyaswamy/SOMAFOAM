/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "emPlasmaFieldsFunctionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineNamedTemplateTypeNameAndDebug(emPlasmaFieldsFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        emPlasmaFieldsFunctionObject,
        dictionary
    );
}

// ************************************************************************* //
