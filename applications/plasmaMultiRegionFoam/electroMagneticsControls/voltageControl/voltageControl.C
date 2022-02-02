/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "voltageControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace emcModels
{
    defineTypeNameAndDebug(voltage, 0);
    addToRunTimeSelectionTable(emcModel, voltage, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emcModels::voltage::voltage
(
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime
)
:
    emcModel(electroMagnetics, mspm, E, runTime),
	mode_(emcModelCoeffs_.lookup("mode")),
    amplitude_("amplitude", dimless, emcModelCoeffs_.lookup("amplitude")),
    frequency_("frequency", dimless, emcModelCoeffs_.lookup("frequency")),
    bias_("bias", dimless, emcModelCoeffs_.lookup("bias"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::emcModels::voltage::~voltage()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::emcModels::voltage::correct(dictionary& voltageDict)
{
	if (mode_ == "continuousFrequencyModulated")
	{
		scalar voltageValue = amplitude_.value()*Foam::cos(2.0*M_PI*frequency_.value()*time_.value()) + bias_.value();

		voltageDict.set("voltage", voltageValue);
	}
	else
	{
        FatalErrorIn("emcModels::voltage::correct(dictionary& voltageDict)")
            << " incorrect mode "
            << exit(FatalError);
	}
}

bool Foam::emcModels::voltage::read(const dictionary& electroMagnetics)
{
    emcModel::read(electroMagnetics);

    emcModelCoeffs_.lookup("amplitude") >> amplitude_.value();
    emcModelCoeffs_.lookup("frequency") >> frequency_.value();
    emcModelCoeffs_.lookup("bias") >> bias_.value();

    return true;
}

// ************************************************************************* //
