/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

InNamespace
    Foam

Description
    Interpolates y values from one curve to another with a different x
    distribution.

SourceFiles
    plasmaInterpolateXY.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaInterpolateXY_H
#define plasmaInterpolateXY_H

#include "scalar.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Field<Type> plasmaInterpolateXY
(
    const scalarField& xNew,
    const scalarField& xOld,
    const Field<Type>& yOld
);


template<class Type>
Type plasmaInterpolateXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "plasmaInterpolateXY.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
