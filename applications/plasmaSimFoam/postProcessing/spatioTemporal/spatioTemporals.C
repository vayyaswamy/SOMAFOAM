#include "spatioTemporal.H"

namespace Foam
{
	defineTypeNameAndDebug(spatioTemporal, 0);
	defineRunTimeSelectionTable(spatioTemporal, dictionary);
};

Foam::spatioTemporal::spatioTemporal
(
	const dictionary& spatioTemporals,
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime,
	const fvMesh& mesh
):
	regIOobject
	(
		IOobject
		(
			"spatioTemporals",
			E.time().constant(),
			E.db()
		)
	),
	spatioTemporalcoeffs_
	(
		spatioTemporals.subDict
		(
			word(spatioTemporals.lookup("spatioTemporal")) + "Coeffs"
		)
	),
	emcModelCoeffs_
	(
		electroMagnetics.subDict
		(
			word(electroMagnetics.lookup("emcModel")) + "Coeffs"
		)
	),
	mspm_(mspm),
	E_(E),
	time_(runTime),
	mesh_(mesh)
	{

	}

Foam::spatioTemporal::~spatioTemporal()
{

}


