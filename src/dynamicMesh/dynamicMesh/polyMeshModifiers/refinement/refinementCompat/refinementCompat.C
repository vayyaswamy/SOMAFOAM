// Functions added to namespace meshTools in foam-extend 4.1
// (src/foam/meshes/meshTools/meshTools.C), see refinementCompat.H

#include "refinementCompat.H"
#include "polyMesh.H"
#include "faceZone.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::meshTools::setFaceInfo
(
    const polyMesh& mesh,
    const label faceI,
    label& patchID,
    label& zoneID,
    label& zoneFlip
)
{
    patchID = -1;

    if (!mesh.isInternalFace(faceI))
    {
        patchID = mesh.boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh.faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID > -1)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}


void Foam::meshTools::extendMarkedCellsAcrossFaces
(
    const polyMesh& mesh,
    boolList& markedCell
)
{
    // Mark all faces for all marked cells
    const label nFaces = mesh.nFaces();
    boolList markedFace(nFaces, false);

    // Get mesh cells
    const cellList& meshCells = mesh.cells();

    // Loop through all cells
    forAll (markedCell, cellI)
    {
        if (markedCell[cellI])
        {
            // This cell is marked, get its faces
            const cell& cFaces = meshCells[cellI];

            forAll (cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    // Snyc the face list across processor boundaries
    syncTools::syncFaceList(mesh, markedFace, orEqOp<bool>(), false);

    // Get necessary mesh data
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    // Internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Face is marked, mark both owner and neighbour
            const label& own = owner[faceI];
            const label& nei = neighbour[faceI];

            // Mark owner and neighbour cells
            markedCell[own] = true;
            markedCell[nei] = true;
        }
    }

    // Boundary faces
    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Face is marked, mark owner
            const label& own = owner[faceI];

            // Mark owner
            markedCell[own] = true;
        }
    }
}


void Foam::meshTools::extendMarkedCellsAcrossPoints
(
    const polyMesh& mesh,
    boolList& markedCell
)
{
    // Mark all points for all marked cells
    const label nPoints = mesh.nPoints();
    boolList markedPoint(nPoints, false);

    // Get cell points
    const labelListList& meshCellPoints = mesh.cellPoints();

    // Loop through all cells
    forAll (markedCell, cellI)
    {
        if (markedCell[cellI])
        {
            // This cell is marked, get its points
            const labelList& cPoints = meshCellPoints[cellI];

            forAll (cPoints, i)
            {
                markedPoint[cPoints[i]] = true;
            }
        }
    }

    // Snyc point list across processor boundaries
    syncTools::syncPointList
    (
        mesh,
        markedPoint,
        orEqOp<bool>(),
        true, // Default value
        true  // Apply separation for parallel cyclics
    );

    // Get point cells
    const labelListList& meshPointCells = mesh.pointCells();

    // Loop through all points
    forAll (markedPoint, pointI)
    {
        if (markedPoint[pointI])
        {
            // This point is marked, mark all of its cells
            const labelList& pCells = meshPointCells[pointI];

            forAll (pCells, i)
            {
                markedCell[pCells[i]] = true;
            }
        }
    }
}


// ************************************************************************* //
