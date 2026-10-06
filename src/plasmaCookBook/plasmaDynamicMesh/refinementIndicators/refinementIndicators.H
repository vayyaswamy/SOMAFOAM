/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Namespace
    Foam::refinementIndicators

Description
    Refinement indicator built from a list of entries, shared by the adaptive
    mesh classes. The indicator of a cell is the maximum over the entries of
    weight*value, with value depending on the entry type:

        relativeGradient  h |grad f|/(|f| + floor): the gradient relative to
                          the value, times the size h of the cell, i.e. the
                          cell size over the local scale length of f
        gradient          jump between neighbouring cells, |f1 - f2|/scale
        magnitude         |f|/scale

    for any volScalarField or volVectorField f of the mesh:

    \verbatim
        indicators
        (
            { type relativeGradient; field N_electron; floor 1e12; }
            { type magnitude;        field E;  scale 1e6; weight 0.5; }
        );
    \endverbatim

SourceFiles
    refinementIndicators.C

\*---------------------------------------------------------------------------*/

#ifndef refinementIndicators_H
#define refinementIndicators_H

#include "fvMesh.H"
#include "dictionary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementIndicators
{
    //- Evaluate the "indicators" entries of dict; indicator is resized to
    //  the number of cells
    void calculate
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalarField& indicator
    );
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
