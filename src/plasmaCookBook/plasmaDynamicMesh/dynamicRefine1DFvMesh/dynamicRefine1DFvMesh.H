/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::dynamicRefine1DFvMesh

Description
    A fvMesh with adaptive refinement along one direction, for 1D meshes
    (a single row of hexahedral cells whose side faces are boundary faces,
    typically on empty patches).

    A cell is refined by splitting it in two along the refinement direction;
    the original cell keeps the low side. Because the side faces belong to
    the cell alone, the mesh stays conformal (no hanging nodes). Two cells
    created by the same split (siblings) can be merged again; merged cells
    take the volume-weighted average of the pair for all volScalarFields and
    volVectorFields.

    Cells are selected from a volScalarField, as in dynamicRefineFvMesh:
    refine where lowerRefineLevel <= field <= upperRefineLevel, merge
    siblings where field < unrefineLevel.

    Example (constant/dynamicMeshDict):
    \verbatim
        dynamicFvMesh   dynamicRefine1DFvMesh;

        dynamicRefine1DFvMeshCoeffs
        {
            direction           (1 0 0);
            refineInterval      10;       // time steps; 0 = never
            field               refinementIndicator;
            lowerRefineLevel    0.5;
            upperRefineLevel    1e30;
            unrefineLevel       0.1;
            nBufferLayers       1;
            maxRefinement       3;
            maxCells            2000;
        }
    \endverbatim

SourceFiles
    dynamicRefine1DFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefine1DFvMesh_H
#define dynamicRefine1DFvMesh_H

#include "dynamicFvMesh.H"
#include "labelIOList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                    Class dynamicRefine1DFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefine1DFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Refinement direction (unit vector)
        vector direction_;

        //- Refinement level of every cell
        labelIOList cellLevel_;

        //- Split that created the cell (-1: original cell). The two cells
        //  created by one split share the value.
        labelIOList cellFamily_;

        //- Per split: level and family of the parent cell (2 entries each)
        labelIOList history_;

        //- Number of mesh changes so far
        label nChanges_;


    // Private Member Functions

        //- The two faces of a cell normal to the refinement direction,
        //  returned as (low side, high side)
        labelPair directionFaces(const label cellI) const;

        //- Split the given cells; returns the number of cells split
        label refine(const labelList& cellsToRefine, boolList& protectedCell);

        //- Merge the sibling cells on both sides of the given faces;
        //  returns the number of faces removed
        label unrefine(const labelList& facesToRemove);

        //- Cells to refine
        labelList selectRefineCells
        (
            const scalarField& vFld,
            const scalar lowerRefineLevel,
            const scalar upperRefineLevel,
            const label maxRefinement,
            const label nBufferLayers,
            const label maxCells
        ) const;

        //- Faces between siblings that may be merged
        labelList selectUnrefineFaces
        (
            const scalarField& vFld,
            const scalar unrefineLevel,
            const boolList& protectedCell
        ) const;

        //- Disallow default bitwise copy construct
        dynamicRefine1DFvMesh(const dynamicRefine1DFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefine1DFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicRefine1DFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefine1DFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicRefine1DFvMesh();


    // Member Functions

        //- Refinement level per cell
        const labelList& cellLevel() const
        {
            return cellLevel_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Write the mesh and the refinement data
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
