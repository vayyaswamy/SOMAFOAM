/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "dynamicRefine1DFvMesh.H"
#include "refinementIndicators.H"
#include "addToRunTimeSelectionTable.H"
#include "directTopoChange.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyAddCell.H"
#include "removeFaces.H"
#include "mapPolyMesh.H"
#include "volFields.H"
#include "HashTable.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefine1DFvMesh, 0);

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicRefine1DFvMesh,
        IOobject
    );


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//- Volume-weighted average over each pair of cells, for all fields of a type
template<class Type>
void storePairAverages
(
    const fvMesh& mesh,
    const List<labelPair>& pairs,
    HashTable<Field<Type>, word>& averages
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    const HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    const scalarField& V = mesh.V().field();

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const Field<Type>& f = iter()->internalField();

        Field<Type> avg(pairs.size());

        forAll(pairs, i)
        {
            const label a = pairs[i].first();
            const label b = pairs[i].second();

            avg[i] = (V[a]*f[a] + V[b]*f[b])/(V[a] + V[b]);
        }

        averages.insert(iter.key(), avg);
    }
}


//- Assign the stored pair averages to the merged cells
template<class Type>
void restorePairAverages
(
    const fvMesh& mesh,
    const labelList& mergedCells,
    const HashTable<Field<Type>, word>& averages
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    forAllConstIter(typename HashTable<Field<Type> >, averages, iter)
    {
        if (mesh.foundObject<GeoField>(iter.key()))
        {
            Field<Type>& f =
                const_cast<GeoField&>
                (
                    mesh.lookupObject<GeoField>(iter.key())
                ).internalField();

            forAll(mergedCells, i)
            {
                f[mergedCells[i]] = iter()[i];
            }
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volScalarField& Foam::dynamicRefine1DFvMesh::calcIndicator
(
    const dictionary& refineDict
)
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
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                *this,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    scalarField& indicator = indicatorPtr_().internalField();

    refinementIndicators::calculate(*this, refineDict, indicator);

    indicatorPtr_().correctBoundaryConditions();

    return indicatorPtr_();
}


Foam::labelPair Foam::dynamicRefine1DFvMesh::directionFaces
(
    const label cellI
) const
{
    const cell& cFaces = cells()[cellI];

    label lowFace = -1;
    label highFace = -1;
    label nFound = 0;

    forAll(cFaces, i)
    {
        const label faceI = cFaces[i];
        const vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);

        if (mag(n & direction_) > 0.5)
        {
            nFound++;

            if (lowFace == -1)
            {
                lowFace = faceI;
            }
            else if
            (
                (faceCentres()[faceI] & direction_)
              < (faceCentres()[lowFace] & direction_)
            )
            {
                highFace = lowFace;
                lowFace = faceI;
            }
            else
            {
                highFace = faceI;
            }
        }
        else if (faceI < nInternalFaces())
        {
            FatalErrorIn("dynamicRefine1DFvMesh::directionFaces(const label)")
                << "Cell " << cellI << " has an internal side face " << faceI
                << ": the mesh is not one-dimensional along " << direction_
                << abort(FatalError);
        }
    }

    if (nFound != 2 || cFaces.size() != 6)
    {
        FatalErrorIn("dynamicRefine1DFvMesh::directionFaces(const label)")
            << "Cell " << cellI << " is not a hexahedron with two faces normal"
            << " to the refinement direction " << direction_
            << " (faces: " << cFaces.size() << ", normal faces: " << nFound
            << ")" << abort(FatalError);
    }

    return labelPair(lowFace, highFace);
}


Foam::label Foam::dynamicRefine1DFvMesh::refine
(
    const labelList& cellsToRefine,
    boolList& protectedCell
)
{
    directTopoChange meshMod(*this);

    const pointField& pts = points();

    forAll(cellsToRefine, i)
    {
        const label cellI = cellsToRefine[i];

        const labelPair dirFaces = directionFaces(cellI);
        const label lowFace = dirFaces.first();
        const label highFace = dirFaces.second();

        const face& fLow = faces()[lowFace];

        labelHashSet lowPoints(2*fLow.size());
        forAll(fLow, fp)
        {
            lowPoints.insert(fLow[fp]);
        }

        // New cell: takes the high side
        const label newCellI = meshMod.setAction
        (
            polyAddCell(-1, -1, -1, cellI, cellZones().whichZone(cellI))
        );

        // Split the side faces; the points added on the edges between the
        // low and the high face are shared by neighbouring side faces
        Map<label> midOfLowPoint(2*fLow.size());
        Map<point> midPosition(2*fLow.size());

        const cell& cFaces = cells()[cellI];

        forAll(cFaces, cf)
        {
            const label faceI = cFaces[cf];

            if (faceI == lowFace || faceI == highFace)
            {
                continue;
            }

            const face& f = faces()[faceI];

            DynamicList<label> lowHalf(f.size());
            DynamicList<label> highHalf(f.size());

            forAll(f, fp)
            {
                const label v = f[fp];
                const label w = f.nextLabel(fp);
                const bool vLow = lowPoints.found(v);

                if (vLow)
                {
                    lowHalf.append(v);
                }
                else
                {
                    highHalf.append(v);
                }

                if (vLow != lowPoints.found(w))
                {
                    const label lowPt = vLow ? v : w;

                    if (!midOfLowPoint.found(lowPt))
                    {
                        const point mid = 0.5*(pts[v] + pts[w]);

                        midOfLowPoint.insert
                        (
                            lowPt,
                            meshMod.setAction
                            (
                                polyAddPoint(mid, lowPt, -1, true)
                            )
                        );
                        midPosition.insert(lowPt, mid);
                    }

                    lowHalf.append(midOfLowPoint[lowPt]);
                    highHalf.append(midOfLowPoint[lowPt]);
                }
            }

            const label patchI = boundaryMesh().whichPatch(faceI);
            const label zoneI = faceZones().whichZone(faceI);
            bool zoneFlip = false;

            if (zoneI >= 0)
            {
                const faceZone& fZone = faceZones()[zoneI];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            lowHalf.shrink();
            highHalf.shrink();

            meshMod.setAction
            (
                polyModifyFace
                (
                    face(lowHalf),
                    faceI,
                    cellI,
                    -1,
                    false,
                    patchI,
                    false,
                    zoneI,
                    zoneFlip
                )
            );

            meshMod.setAction
            (
                polyAddFace
                (
                    face(highHalf),
                    newCellI,
                    -1,
                    -1,
                    -1,
                    faceI,
                    false,
                    patchI,
                    zoneI,
                    zoneFlip
                )
            );
        }

        if (midOfLowPoint.size() != fLow.size())
        {
            FatalErrorIn("dynamicRefine1DFvMesh::refine(const labelList&)")
                << "Cell " << cellI << ": found " << midOfLowPoint.size()
                << " edges along the refinement direction, expected "
                << fLow.size() << abort(FatalError);
        }

        // New internal face between the two halves, pointing from the
        // original (low) cell to the new (high) cell
        {
            face midFace(fLow.size());
            pointField midPts(fLow.size());

            forAll(fLow, fp)
            {
                midFace[fp] = midOfLowPoint[fLow[fp]];
                midPts[fp] = midPosition[fLow[fp]];
            }

            vector n = vector::zero;
            const point c = average(midPts);

            forAll(midPts, fp)
            {
                n += (midPts[fp] - c) ^ (midPts[(fp + 1) % midPts.size()] - c);
            }

            if ((n & direction_) < 0)
            {
                midFace = midFace.reverseFace();
            }

            // Take the flux of an existing internal face as first guess
            label masterFace = -1;

            if (lowFace < nInternalFaces())
            {
                masterFace = lowFace;
            }
            else if (highFace < nInternalFaces())
            {
                masterFace = highFace;
            }

            meshMod.setAction
            (
                polyAddFace
                (
                    midFace,
                    cellI,
                    newCellI,
                    (masterFace == -1 ? fLow[0] : -1),
                    -1,
                    masterFace,
                    false,
                    -1,
                    -1,
                    false
                )
            );
        }

        // The high face now belongs to the new cell
        {
            const label zoneI = faceZones().whichZone(highFace);
            bool zoneFlip = false;

            if (zoneI >= 0)
            {
                const faceZone& fZone = faceZones()[zoneI];
                zoneFlip = fZone.flipMap()[fZone.whichFace(highFace)];
            }

            if (highFace >= nInternalFaces())
            {
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces()[highFace],
                        highFace,
                        newCellI,
                        -1,
                        false,
                        boundaryMesh().whichPatch(highFace),
                        false,
                        zoneI,
                        zoneFlip
                    )
                );
            }
            else if (faceOwner()[highFace] == cellI)
            {
                // The new cell has the highest label, so it becomes the
                // neighbour and the face is turned around
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces()[highFace].reverseFace(),
                        highFace,
                        faceNeighbour()[highFace],
                        newCellI,
                        true,
                        -1,
                        false,
                        zoneI,
                        (zoneI >= 0 ? !zoneFlip : false)
                    )
                );
            }
            else
            {
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces()[highFace],
                        highFace,
                        faceOwner()[highFace],
                        newCellI,
                        false,
                        -1,
                        false,
                        zoneI,
                        zoneFlip
                    )
                );
            }
        }

        // Both halves are siblings of a new family, one level up
        history_.setSize(history_.size() + 2);
        history_[history_.size() - 2] = cellLevel_[cellI];
        history_[history_.size() - 1] = cellFamily_[cellI];

        cellLevel_[cellI]++;
        cellFamily_[cellI] = history_.size()/2 - 1;

        protectedCell[cellI] = true;
    }

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    updateMesh(map());

    // The new cell inherits level, family and protection from its master
    const labelList& cellMap = map().cellMap();

    labelList newLevel(nCells());
    labelList newFamily(nCells());
    boolList newProtected(nCells());

    forAll(cellMap, cellI)
    {
        newLevel[cellI] = cellLevel_[cellMap[cellI]];
        newFamily[cellI] = cellFamily_[cellMap[cellI]];
        newProtected[cellI] = protectedCell[cellMap[cellI]];
    }

    cellLevel_ = newLevel;
    cellFamily_ = newFamily;
    protectedCell.transfer(newProtected);

    nChanges_++;

    return cellsToRefine.size();
}


Foam::label Foam::dynamicRefine1DFvMesh::unrefine
(
    const labelList& facesToRemove
)
{
    // Merge the boundary faces of the merged cells wherever possible
    removeFaces faceRemover(*this, GREAT);

    labelList cellRegion;
    labelList cellRegionMaster;
    labelList allFacesToRemove;

    faceRemover.compatibleRemoves
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        allFacesToRemove
    );

    if (allFacesToRemove.size() != facesToRemove.size())
    {
        FatalErrorIn("dynamicRefine1DFvMesh::unrefine(const labelList&)")
            << "Removing " << facesToRemove.size() << " faces between sibling"
            << " cells requires removing " << allFacesToRemove.size()
            << " faces: the mesh is not one-dimensional"
            << abort(FatalError);
    }

    directTopoChange meshMod(*this);

    // Dummy pointRegionMaster, ignored for cell merging
    labelList pointRegionMaster(cellRegionMaster.size(), label(-1));

    faceRemover.setRefinement
    (
        allFacesToRemove,
        cellRegion,
        pointRegionMaster,
        cellRegionMaster,
        meshMod
    );

    // The surviving cell of each pair (lowest label) and its sibling
    List<labelPair> pairs(facesToRemove.size());

    forAll(facesToRemove, i)
    {
        const label own = faceOwner()[facesToRemove[i]];
        const label nei = faceNeighbour()[facesToRemove[i]];

        pairs[i] = labelPair(min(own, nei), max(own, nei));
    }

    HashTable<scalarField, word> scalarAverages;
    HashTable<vectorField, word> vectorAverages;

    storePairAverages<scalar>(*this, pairs, scalarAverages);
    storePairAverages<vector>(*this, pairs, vectorAverages);

    // Fields averaged with an additional weight, e.g. Te with the electron
    // density; the old-time fields are treated alike
    {
        const scalarField& V = this->V().field();

        forAll(weightedAverages_, wI)
        {
            for (label old = 0; old < 2; old++)
            {
                const word fName
                (
                    weightedAverages_[wI].first() + (old ? "_0" : "")
                );
                const word wName
                (
                    weightedAverages_[wI].second() + (old ? "_0" : "")
                );

                if
                (
                    !foundObject<volScalarField>(fName)
                 || !foundObject<volScalarField>(wName)
                )
                {
                    if (!old)
                    {
                        FatalErrorIn
                        (
                            "dynamicRefine1DFvMesh::unrefine(const labelList&)"
                        )   << "weightedAverages: field " << fName
                            << " or weight " << wName << " not found"
                            << abort(FatalError);
                    }
                    continue;
                }

                const scalarField& f =
                    lookupObject<volScalarField>(fName).internalField();
                const scalarField& w =
                    lookupObject<volScalarField>(wName).internalField();

                scalarField& avg = scalarAverages[fName];

                forAll(pairs, i)
                {
                    const label a = pairs[i].first();
                    const label b = pairs[i].second();
                    const scalar wa = V[a]*mag(w[a]);
                    const scalar wb = V[b]*mag(w[b]);

                    if (wa + wb > VSMALL)
                    {
                        avg[i] = (wa*f[a] + wb*f[b])/(wa + wb);
                    }
                }
            }
        }
    }

    // The merged cell gets the level and family of the parent
    forAll(pairs, i)
    {
        const label family = cellFamily_[pairs[i].first()];
        const label parentLevel = history_[2*family];
        const label parentFamily = history_[2*family + 1];

        cellLevel_[pairs[i].first()] = parentLevel;
        cellFamily_[pairs[i].first()] = parentFamily;
        cellLevel_[pairs[i].second()] = parentLevel;
        cellFamily_[pairs[i].second()] = parentFamily;
    }

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    updateMesh(map());

    const labelList& cellMap = map().cellMap();
    const labelList& reverseCellMap = map().reverseCellMap();

    labelList newLevel(nCells());
    labelList newFamily(nCells());

    forAll(cellMap, cellI)
    {
        newLevel[cellI] = cellLevel_[cellMap[cellI]];
        newFamily[cellI] = cellFamily_[cellMap[cellI]];
    }

    cellLevel_ = newLevel;
    cellFamily_ = newFamily;

    labelList mergedCells(pairs.size());

    forAll(pairs, i)
    {
        mergedCells[i] = reverseCellMap[pairs[i].first()];

        if (mergedCells[i] < 0)
        {
            FatalErrorIn("dynamicRefine1DFvMesh::unrefine(const labelList&)")
                << "Cell " << pairs[i].first() << " did not survive the merge"
                << abort(FatalError);
        }
    }

    restorePairAverages<scalar>(*this, mergedCells, scalarAverages);
    restorePairAverages<vector>(*this, mergedCells, vectorAverages);

    nChanges_++;

    return facesToRemove.size();
}


Foam::labelList Foam::dynamicRefine1DFvMesh::selectRefineCells
(
    const scalarField& vFld,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const label maxRefinement,
    const label nBufferLayers,
    const label maxCells
) const
{
    boolList marked(nCells(), false);

    forAll(vFld, cellI)
    {
        if
        (
            vFld[cellI] >= lowerRefineLevel
         && vFld[cellI] <= upperRefineLevel
         && cellLevel_[cellI] < maxRefinement
        )
        {
            marked[cellI] = true;
        }
    }

    const labelListList& cc = cellCells();

    // Buffer layers around the marked cells
    for (label layer = 0; layer < nBufferLayers; layer++)
    {
        boolList extended(marked);

        forAll(marked, cellI)
        {
            if (marked[cellI])
            {
                forAll(cc[cellI], i)
                {
                    const label nbr = cc[cellI][i];

                    if (cellLevel_[nbr] < maxRefinement)
                    {
                        extended[nbr] = true;
                    }
                }
            }
        }

        marked.transfer(extended);
    }

    // Keep the level difference between neighbours at one: a coarser
    // neighbour of a marked cell is refined as well
    bool changed = true;

    while (changed)
    {
        changed = false;

        forAll(marked, cellI)
        {
            if (marked[cellI])
            {
                forAll(cc[cellI], i)
                {
                    const label nbr = cc[cellI][i];

                    if (!marked[nbr] && cellLevel_[nbr] < cellLevel_[cellI])
                    {
                        marked[nbr] = true;
                        changed = true;
                    }
                }
            }
        }
    }

    const label nAllowed = max(maxCells - nCells(), 0);

    DynamicList<label> cellsToRefine(nCells());

    forAll(marked, cellI)
    {
        if (marked[cellI] && cellsToRefine.size() < nAllowed)
        {
            cellsToRefine.append(cellI);
        }
    }

    cellsToRefine.shrink();

    return labelList(cellsToRefine);
}


Foam::labelList Foam::dynamicRefine1DFvMesh::selectUnrefineFaces
(
    const scalarField& vFld,
    const scalar unrefineLevel,
    const boolList& protectedCell
) const
{
    const labelListList& cc = cellCells();

    DynamicList<label> facesToRemove(nInternalFaces());

    for (label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        const label own = faceOwner()[faceI];
        const label nei = faceNeighbour()[faceI];

        if
        (
            cellFamily_[own] < 0
         || cellFamily_[own] != cellFamily_[nei]
         || cellLevel_[own] != cellLevel_[nei]
         || protectedCell[own]
         || protectedCell[nei]
         || vFld[own] >= unrefineLevel
         || vFld[nei] >= unrefineLevel
        )
        {
            continue;
        }

        // Merging must not leave a neighbour two levels finer
        bool finerNeighbour = false;

        forAll(cc[own], i)
        {
            if (cellLevel_[cc[own][i]] > cellLevel_[own])
            {
                finerNeighbour = true;
            }
        }

        forAll(cc[nei], i)
        {
            if (cellLevel_[cc[nei][i]] > cellLevel_[nei])
            {
                finerNeighbour = true;
            }
        }

        if (!finerNeighbour)
        {
            facesToRemove.append(faceI);
        }
    }

    facesToRemove.shrink();

    return labelList(facesToRemove);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefine1DFvMesh::dynamicRefine1DFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    direction_(1, 0, 0),
    cellLevel_
    (
        IOobject
        (
            "cellLevel1D",
            time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(nCells(), 0)
    ),
    cellFamily_
    (
        IOobject
        (
            "cellFamily1D",
            time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(nCells(), -1)
    ),
    history_
    (
        IOobject
        (
            "refinementHistory1D",
            time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(0)
    ),
    nChanges_(0),
    indicatorPtr_(),
    weightedAverages_()
{
    const dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    );

    direction_ = refineDict.lookupOrDefault<vector>("direction", vector(1, 0, 0));
    direction_ /= mag(direction_);

    if (cellLevel_.size() != nCells() || cellFamily_.size() != nCells())
    {
        FatalErrorIn("dynamicRefine1DFvMesh::dynamicRefine1DFvMesh(const IOobject&)")
            << "Refinement data read for " << cellLevel_.size()
            << " cells, but the mesh has " << nCells() << " cells"
            << abort(FatalError);
    }

    // Check that the mesh is one-dimensional along the direction
    forAll(cells(), cellI)
    {
        directionFaces(cellI);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefine1DFvMesh::~dynamicRefine1DFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefine1DFvMesh::update()
{
    // Re-read the dictionary so the controls can be changed during the run
    const dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    );

    const label refineInterval = readLabel(refineDict.lookup("refineInterval"));

    bool hasChanged = false;

    if
    (
        refineInterval > 0
     && time().timeIndex() > 0
     && time().timeIndex() % refineInterval == 0
    )
    {
        const label maxCells = readLabel(refineDict.lookup("maxCells"));
        const label maxRefinement =
            readLabel(refineDict.lookup("maxRefinement"));
        const scalar lowerRefineLevel =
            readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel =
            readScalar(refineDict.lookup("upperRefineLevel"));
        const scalar unrefineLevel =
            readScalar(refineDict.lookup("unrefineLevel"));
        const label nBufferLayers =
            readLabel(refineDict.lookup("nBufferLayers"));

        weightedAverages_ = refineDict.lookupOrDefault<List<Pair<word> > >
        (
            "weightedAverages",
            List<Pair<word> >(0)
        );

        // Indicator: built here from the indicators entries, or a field
        // provided by the solver
        const volScalarField& vFld =
        (
            refineDict.found("indicators")
          ? calcIndicator(refineDict)
          : lookupObject<volScalarField>(word(refineDict.lookup("field")))
        );

        // Cells refined in this update are not merged again straight away
        boolList protectedCell(nCells(), false);

        const label nOldCells = nCells();

        const labelList cellsToRefine
        (
            selectRefineCells
            (
                vFld.internalField(),
                lowerRefineLevel,
                upperRefineLevel,
                maxRefinement,
                nBufferLayers,
                maxCells
            )
        );

        if (cellsToRefine.size())
        {
            refine(cellsToRefine, protectedCell);

            Info<< "dynamicRefine1DFvMesh: refined " << cellsToRefine.size()
                << " cells, " << nOldCells << " -> " << nCells() << " cells"
                << endl;

            hasChanged = true;
        }

        // Cells in the buffer layers around cells that call for refinement
        // stay refined, otherwise they would be merged and split again at
        // every update
        {
            const scalarField& v = vFld.internalField();
            const labelListList& cc = cellCells();

            boolList keep(nCells(), false);

            forAll(v, cellI)
            {
                if (v[cellI] >= lowerRefineLevel && v[cellI] <= upperRefineLevel)
                {
                    keep[cellI] = true;
                }
            }

            for (label layer = 0; layer < nBufferLayers; layer++)
            {
                boolList extended(keep);

                forAll(keep, cellI)
                {
                    if (keep[cellI])
                    {
                        forAll(cc[cellI], i)
                        {
                            extended[cc[cellI][i]] = true;
                        }
                    }
                }

                keep.transfer(extended);
            }

            forAll(keep, cellI)
            {
                if (keep[cellI])
                {
                    protectedCell[cellI] = true;
                }
            }
        }

        const labelList facesToRemove
        (
            selectUnrefineFaces
            (
                vFld.internalField(),
                unrefineLevel,
                protectedCell
            )
        );

        if (facesToRemove.size())
        {
            const label nBefore = nCells();

            unrefine(facesToRemove);

            Info<< "dynamicRefine1DFvMesh: merged " << facesToRemove.size()
                << " cell pairs, " << nBefore << " -> " << nCells()
                << " cells" << endl;

            hasChanged = true;
        }
    }

    changing(hasChanged);

    return hasChanged;
}


bool Foam::dynamicRefine1DFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    bool ok = dynamicFvMesh::writeObject(fmt, ver, cmp);

    // The refinement data go with the fields, so a restart picks them up
    if (nChanges_ > 0)
    {
        const_cast<labelIOList&>(cellLevel_).instance() = time().timeName();
        const_cast<labelIOList&>(cellFamily_).instance() = time().timeName();
        const_cast<labelIOList&>(history_).instance() = time().timeName();

        ok = ok
          && cellLevel_.writeObject(fmt, ver, cmp)
          && cellFamily_.writeObject(fmt, ver, cmp)
          && history_.writeObject(fmt, ver, cmp);
    }

    return ok;
}


// ************************************************************************* //
