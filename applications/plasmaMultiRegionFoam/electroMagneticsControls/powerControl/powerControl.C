/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "powerControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace emcModels
{
    defineTypeNameAndDebug(power, 0);
    addToRunTimeSelectionTable(emcModel, power, dictionary);
};
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emcModels::power::power
(
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime
)
:
    emcModel(electroMagnetics, mspm, E, runTime),
	mode_(emcModelCoeffs_.lookup("mode")),
    initialAmplitude_(readScalar(emcModelCoeffs_.lookup("initialAmplitude"))),
    frequency_(readScalar(emcModelCoeffs_.lookup("frequency"))),
    bias_(readScalar(emcModelCoeffs_.lookup("bias"))),
    power_(readScalar(emcModelCoeffs_.lookup("power"))),
    nCycles_(readInt(emcModelCoeffs_.lookup("controlFrequency"))),
    dampingFactor_(readScalar(emcModelCoeffs_.lookup("dampingFactor"))),
    mf_(readScalar(emcModelCoeffs_.lookup("geometricFactor"))),
	timeCounter_(0.0),
	timeCount_(0.0),
	curTimeIndex_(time_.timeIndex()),
	powerLogFilePtr_(NULL),
    amplitude_(0.0),
    powerSum_(0.0),
    meshV_
    (
        IOobject
        (
            "meshV",
            E_.mesh().time().timeName(),
            E_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        E_.mesh(),
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0)
    )
{
	meshV_.internalField() = E_.mesh().V();

    if (Pstream::master())
    {
        powerLogFilePtr_ =
            new OFstream
            (fileName("power_voltage_logfile"));
        OFstream& powerLogFile = *powerLogFilePtr_;
        int width = 20;
        powerLogFile << "time";
        powerLogFile.width(width);
        powerLogFile << "power";
        powerLogFile.width(width);
        powerLogFile << "voltage";
		powerLogFile << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::emcModels::power::~power()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::emcModels::power::powerSumMesh() const
{
    const objectRegistry& db = E_.db();
    const volVectorField& tddtE = db.lookupObject<volVectorField>("ddtE");

	volScalarField tpowerSumMesh = meshV_*((plasmaConstants::eCharge*mspm_.netChargeFlux() + plasmaConstants::epsilon0*tddtE) & E_);

    return gSum(tpowerSumMesh);
}

void Foam::emcModels::power::correct(dictionary& voltageDict)
{
	if (mode_ == "continuousFrequencyModulated")
	{
		const scalar& tpower = powerSumMesh();

		powerSum_ += tpower;

		curTimeIndex_ = time_.timeIndex();

		if(curTimeIndex_ == 1)
		{
			amplitude_ = initialAmplitude_;
		}

		timeCount_ = nCycles_/frequency_/time_.deltaT().value();

		if(curTimeIndex_ - timeCounter_ >= timeCount_)
		{
			timeCounter_ =  curTimeIndex_;

			scalar powerSumAve_ = powerSum_*mf_/timeCount_;

			powerSum_ = 0.0;

			scalar amplitudeold_(amplitude_);

			amplitude_ = amplitudeold_*(1.0-dampingFactor_*((powerSumAve_/power_)-1.0));

			if((amplitude_/amplitudeold_) >= 1.1)
			{
				amplitude_ = 1.1*amplitudeold_;
			}
			else if((amplitude_/amplitudeold_) <= 0.9)
			{
				amplitude_ = 0.9*amplitudeold_;
			}

			if (Pstream::master())
			{
			   OFstream& powerLogFile = *powerLogFilePtr_;
			   int width = 20;
			   powerLogFile << time_.value();
			   powerLogFile.width(width);
			   powerLogFile << powerSumAve_;
			   powerLogFile.width(width);
			   powerLogFile << amplitude_;
			   powerLogFile << endl;
			}
		}

		scalar voltageValue = amplitude_*Foam::cos(2.0*M_PI*frequency_*time_.value()) + bias_;

		voltageDict.set("voltage", voltageValue);
	}
	else
	{
        FatalErrorIn("emcModels::power::correct(dictionary& voltageDict)")
            << " incorrect mode "
            << exit(FatalError);
	}
}


bool Foam::emcModels::power::read(const dictionary& electroMagnetics)
{
    emcModel::read(electroMagnetics);

    emcModelCoeffs_.lookup("initialAmplitude") >> initialAmplitude_;
    emcModelCoeffs_.lookup("frequency") >> frequency_;
    emcModelCoeffs_.lookup("bias") >> bias_;
    emcModelCoeffs_.lookup("power") >> power_;
    emcModelCoeffs_.lookup("controlFrequency") >> nCycles_;
    emcModelCoeffs_.lookup("dampingFactor") >> dampingFactor_;
    emcModelCoeffs_.lookup("geometricFactor") >> mf_;

    return true;
}


// ************************************************************************* //
