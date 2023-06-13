#include "spatioTemporal.H"

//  
Foam::autoPtr<Foam::spatioTemporal> Foam::spatioTemporal::New
(
	const dictionary& spatioTemporals,
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime,
	const fvMesh& mesh
)
{
	// main dictionary to define subDictionary
	word spatioTemporalTypeName =  spatioTemporals.lookup("spatioTemporal");

	Info << "Initializing 1D spatio-temporal tool " << endl;

	// Type of spatio-temporal tool 
	if (spatioTemporalTypeName == "noneD")
	{
		Info << "1D spatio-temporal tool was not initialized " << nl << endl;
	}

	dictionaryConstructorTable::iterator cstrIter =
		dictionaryConstructorTablePtr_->find(spatioTemporalTypeName);

	// none was defined 
	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"spatioTemporal::New"
		)	<< "Unknown spatioTemporal type "
			<< spatioTemporalTypeName << endl << endl
			<< "Valid spatioTemporal are: " << endl
			<< dictionaryConstructorTablePtr_->toc()
			<< exit(FatalError);
	}
	
	return autoPtr<spatioTemporal>
		(cstrIter()(spatioTemporals, electroMagnetics, mspm, E, runTime, mesh));
}