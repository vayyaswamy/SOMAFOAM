#include "spatioNone.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
	namespace spatioTemporals
	{
		defineTypeNameAndDebug(spatioNone, 0);
		addToRunTimeSelectionTable(spatioTemporal, spatioNone, dictionary);
	};
};

Foam::spatioTemporals::spatioNone::spatioNone
(
	const dictionary& spatioTemporals,
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime,
	const fvMesh& mesh
)
:
	spatioTemporal(spatioTemporals,electroMagnetics,mspm,E,runTime,mesh)
	{

	}

Foam::spatioTemporals::spatioNone::~spatioNone()
{

}

void Foam::spatioTemporals::spatioNone::correct()
{
	
}