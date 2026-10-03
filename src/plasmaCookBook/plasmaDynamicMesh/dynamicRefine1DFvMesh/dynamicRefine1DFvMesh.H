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

    Cells are selected from a refinement indicator: refine where
    lowerRefineLevel <= indicator <= upperRefineLevel, merge siblings where
    indicator < unrefineLevel. The indicator is either an existing
    volScalarField (keyword field, as in dynamicRefineFvMesh) or is built
    from a list of entries (keyword indicators) as the maximum over the
    entries of weight*value, with value depending on the entry type:

        relativeGradient  jump between neighbouring cells relative to their
                          mean, |f1 - f2|/(0.5(|f1| + |f2|) + floor);
                          the relative gradient times the cell size
        gradient          |f1 - f2|/scale
        magnitude         |f|/scale

    The fields may be volScalarFields or volVectorFields. The indicator is
    stored as the volScalarField refinementIndicator and written with the
    other fields.

    Fields listed in weightedAverages as (field weightField) are averaged
    with weight volume*weightField when cells merge, e.g. (Te N_electron)
    to conserve the electron thermal energy.

    Example (constant/dynamicMeshDict):
    \verbatim
        dynamicFvMesh   dynamicRefine1DFvMesh;

        dynamicRefine1DFvMeshCoeffs
        {
            direction           (1 0 0);
            refineInterval      10;       // time steps; 0 = never
            indicators
            (
                { type relativeGradient; field N_electron; floor 1e14; }
                { type relativeGradient; field N_Arp1; weight 0.5; }
                { type magnitude;        field E;  scale 1e6; }
            );
            weightedAverages    ((Te N_electron));
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
#include "volFieldsFwd.H"
#include "Pair.H"
#include "autoPtr.H"

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

        //- Refinement indicator built from the indicators entries
        mutable autoPtr<volScalarField> indicatorPtr_;

        //- (field weightField) pairs averaged with weight on merging
        List<Pair<word> > weightedAverages_;


    // Private Member Functions

        //- The two faces of a cell normal to the refinement direction,
        //  returned as (low side, high side)
        labelPair directionFaces(const label cellI) const;

        //- Split the given cells; returns the number of cells split
        label refine(const labelList& cellsToRefine, boolList& protectedCell);

        //- Merge the sibling cells on both sides of the given faces;
        //  returns the number of faces removed
        label unrefine(const labelList& facesToRemove);

        //- Build the refinement indicator from the indicators entries
        const volScalarField& calcIndicator(const dictionary& refineDict);

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
