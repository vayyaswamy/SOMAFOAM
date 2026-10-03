/*---------------------------------------------------------------------------*\
Description
    Compatibility layer for the adaptive refinement classes ported from
    foam-extend 4.1: functions that 4.1 added to namespace meshTools and the
    ...InFunction message macros. Included ahead of the ported sources through
    the -include flag in Make/options.
\*---------------------------------------------------------------------------*/

#ifndef refinementCompat_H
#define refinementCompat_H

#include "meshTools.H"
#include "boolList.H"
#include "error.H"

#ifndef FUNCTION_NAME
#define FUNCTION_NAME __PRETTY_FUNCTION__
#endif
#ifndef FatalErrorInFunction
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)
#endif
#ifndef FatalIOErrorInFunction
#define FatalIOErrorInFunction(ios) FatalIOErrorIn(FUNCTION_NAME, ios)
#endif
#ifndef WarningInFunction
#define WarningInFunction WarningIn(FUNCTION_NAME)
#endif
#ifndef InfoInFunction
#define InfoInFunction InfoIn(FUNCTION_NAME)
#endif

namespace Foam
{

class polyMesh;

namespace meshTools
{
    //- Set face information: patch, zone and zone flip for a face
    void setFaceInfo
    (
        const polyMesh& mesh,
        const label faceI,
        label& patchID,
        label& zoneID,
        label& zoneFlip
    );

    //- Extend marked cells across faces
    void extendMarkedCellsAcrossFaces
    (
        const polyMesh& mesh,
        boolList& markedCell
    );

    //- Extend marked cells across points
    void extendMarkedCellsAcrossPoints
    (
        const polyMesh& mesh,
        boolList& markedCell
    );
}

}

#endif
