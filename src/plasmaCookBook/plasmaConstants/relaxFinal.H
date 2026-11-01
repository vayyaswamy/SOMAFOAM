/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Global
    Foam::relaxFinal

Description
    Field relaxation that uses the "<field>Final" relaxation factor of
    fvSolution on the final PIMPLE iteration of a time step, and the factor
    of "<field>" otherwise. GeometricField::relax() of foam-extend 4.0 always
    uses the latter, so that with

        fields { ".*" 0.8; ".*Final" 1.0; }

    and a single PIMPLE iteration every time step would take only 80 % of
    its change.

    finalIteration = true selects the final factor for a solve that is made
    once per time step outside the PIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef relaxFinal_H
#define relaxFinal_H

#include "GeometricField.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, template<class> class PatchField, class GeoMesh>
void relaxFinal
(
    GeometricField<Type, PatchField, GeoMesh>& field,
    const bool finalIteration = false
)
{
    const solution& sol = field.mesh().solutionDict();

    const word finalName(field.name() + "Final");

    const bool useFinal =
        (finalIteration || sol.lookupOrDefault<bool>("finalIteration", false))
     && sol.relaxField(finalName);

    if (!useFinal && !sol.relaxField(field.name()))
    {
        return;
    }

    const scalar alpha =
        sol.fieldRelaxationFactor(useFinal ? finalName : field.name());

    // A factor of 1 leaves the field as it is: the previous iteration is
    // then not needed (somaFoam does not store it for every field)
    if (alpha != 1)
    {
        field.relax(alpha);
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
