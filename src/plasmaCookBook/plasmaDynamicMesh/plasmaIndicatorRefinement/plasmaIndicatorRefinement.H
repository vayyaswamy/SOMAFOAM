/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::plasmaIndicatorRefinement

Description
    Refinement selection for dynamicPolyRefinementFvMesh (2D and 3D) based on
    the refinement indicator of refinementIndicators: cells are candidates
    for refinement where lowerRefineLevel <= indicator <= upperRefineLevel,
    and a split point is a candidate for unrefinement when the indicator of
    all cells around it is below unrefineLevel.

    The indicator is stored as the volScalarField refinementIndicator.

    Example (constant/dynamicMeshDict):
    \verbatim
        dynamicFvMesh   dynamicPolyRefinementFvMesh;

        dynamicPolyRefinementFvMeshCoeffs
        {
            refineInterval      50;
            unrefineInterval    50;
            separateUpdates     false;

            active              yes;
            maxCells            200000;
            maxRefinementLevel  2;
            nRefinementBufferLayers   2;
            nUnrefinementBufferLayers 4;
            edgeBasedConsistency yes;

            refinementSelection
            {
                type            plasmaIndicatorRefinement;
                indicators
                (
                    { type relativeGradient; field N_electron; floor 1e12; }
                    { type relativeGradient; field N_Arp1;     floor 1e12; }
                );
                lowerRefineLevel    0.2;
                upperRefineLevel    1e30;   // optional
                unrefineLevel       0.07;
            }
        }
    \endverbatim

SourceFiles
    plasmaIndicatorRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaIndicatorRefinement_H
#define plasmaIndicatorRefinement_H

#include "refinementSelection.H"
#include "volFieldsFwd.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class plasmaIndicatorRefinement Declaration
\*---------------------------------------------------------------------------*/

class plasmaIndicatorRefinement
:
    public refinementSelection
{
    // Private data

        //- Refine where the indicator is between these levels
        scalar lowerRefineLevel_;
        scalar upperRefineLevel_;

        //- Unrefine where the indicator is below this level
        scalar unrefineLevel_;

        //- The indicator, kept as a field for output
        mutable autoPtr<volScalarField> indicatorPtr_;


    // Private Member Functions

        //- Evaluate the indicator on the current mesh
        const volScalarField& indicator() const;

        //- Disallow default bitwise copy construct
        plasmaIndicatorRefinement(const plasmaIndicatorRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const plasmaIndicatorRefinement&);


public:

    //- Runtime type information
    TypeName("plasmaIndicatorRefinement");


    // Constructors

        //- Construct from components
        plasmaIndicatorRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~plasmaIndicatorRefinement();


    // Member Functions

        //- Return transferable list of cells to refine
        virtual Xfer<labelList> refinementCellCandidates() const;

        //- Return transferable list of split points to unrefine
        virtual Xfer<labelList> unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
