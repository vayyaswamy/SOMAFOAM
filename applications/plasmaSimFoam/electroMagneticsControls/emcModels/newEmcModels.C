/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "emcModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::emcModel> Foam::emcModel::New
(
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime
)
{
    word emcModelTypeName = electroMagnetics.lookup("emcModel");

    Info<< "Selecting electromagnetic control " << emcModelTypeName << endl;

	if (emcModelTypeName == "none")
	{
		    Info<< "electromagnetic control model " << emcModelTypeName
        << " not implemented." << nl << endl;
	}

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(emcModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "emcModel::New"
        )   << "Unknown emcModel type "
            << emcModelTypeName << endl << endl
            << "Valid emcModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<emcModel>
        (cstrIter()(electroMagnetics, mspm, E, runTime));
}


// ************************************************************************* //
