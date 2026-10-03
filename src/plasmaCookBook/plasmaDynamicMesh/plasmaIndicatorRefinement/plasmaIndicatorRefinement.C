/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaIndicatorRefinement.H"
#include "refinementIndicators.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plasmaIndicatorRefinement, 0);

    addToRunTimeSelectionTable
    (
        refinementSelection,
        plasmaIndicatorRefinement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volScalarField& Foam::plasmaIndicatorRefinement::indicator() const
{
    if (!indicatorPtr_.valid())
    {
        indicatorPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "refinementIndicator",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    refinementIndicators::calculate
    (
        mesh(),
        coeffDict(),
        indicatorPtr_().internalField()
    );

    indicatorPtr_().correctBoundaryConditions();

    return indicatorPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaIndicatorRefinement::plasmaIndicatorRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    lowerRefineLevel_(readScalar(coeffDict().lookup("lowerRefineLevel"))),
    upperRefineLevel_
    (
        coeffDict().lookupOrDefault<scalar>("upperRefineLevel", GREAT)
    ),
    unrefineLevel_(readScalar(coeffDict().lookup("unrefineLevel"))),
    indicatorPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::plasmaIndicatorRefinement::~plasmaIndicatorRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::plasmaIndicatorRefinement::refinementCellCandidates() const
{
    const scalarField& ind = indicator().internalField();

    DynamicList<label> candidates(Foam::max(100, mesh().nCells()/5));

    forAll(ind, cellI)
    {
        if (ind[cellI] >= lowerRefineLevel_ && ind[cellI] <= upperRefineLevel_)
        {
            candidates.append(cellI);
        }
    }

    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(candidates.size(), sumOp<label>())
        << " cells as refinement candidates." << endl;

    return candidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::plasmaIndicatorRefinement::unrefinementPointCandidates() const
{
    const scalarField& ind = indicator().internalField();

    const labelListList& pointCells = mesh().pointCells();

    DynamicList<label> candidates(Foam::max(100, mesh().nPoints()/10));

    // A point qualifies when all the cells around it are below the level
    forAll(pointCells, pointI)
    {
        const labelList& pCells = pointCells[pointI];

        bool low = true;

        forAll(pCells, i)
        {
            if (ind[pCells[i]] >= unrefineLevel_)
            {
                low = false;
                break;
            }
        }

        if (low)
        {
            candidates.append(pointI);
        }
    }

    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(candidates.size(), sumOp<label>())
        << " points as unrefinement candidates." << endl;

    return candidates.xfer();
}


// ************************************************************************* //
