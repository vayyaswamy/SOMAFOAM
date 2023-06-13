/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/
#include "currentControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
	namespace emcModels
	{
		defineTypeNameAndDebug(current, 0);
		addToRunTimeSelectionTable(emcModel, current, dictionary);
	};
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::emcModels::current::current
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
rms_(readScalar(emcModelCoeffs_.lookup("rms"))),
dampingFactor_(readScalar(emcModelCoeffs_.lookup("dampingFactor"))),
timeCount_(0.0),
curTimeIndex_(time_.timeIndex()),
LogFilePtr_(NULL),
amplitude_(0.0),
rmsCurrent_(0.0),
ncells_(0.0),
currentSum_(0.0),
operation_(emcModelCoeffs_.lookup("waveform")),
waveform_(emcModelCoeffs_.lookup("operation")),
dutyCycle_(readScalar(emcModelCoeffs_.lookup("dutyCycle"))),
w_(readScalar(emcModelCoeffs_.lookup("naturalFrequency"))),
e_(readScalar(emcModelCoeffs_.lookup("dampingRatio"))),
powerSum_(0.0),
tolerance_(readScalar(emcModelCoeffs_.lookup("tolerance"))),
tpowerOld_(0.0),
tcurrentDensityOld_(0.0),
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

	forAll(E_,celli)
	{
		ncells_ = ncells_ + 1;
	}

    if (Pstream::master())
	{
    	LogFilePtr_ = new OFstream(fileName("current_voltage_logfile"));
    	OFstream& resistanceLogFile = *LogFilePtr_;
    	resistanceLogFile << "time" << tab;
    	resistanceLogFile << "input" << tab;
    	resistanceLogFile << "RMScurrent" << tab;
    	resistanceLogFile << "AVGpower" << endl;
	}
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //
Foam::emcModels::current::~current()
{

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::scalar Foam::emcModels::current::currentDensitySum() const
{
	const objectRegistry& db = E_.db();
	const volVectorField& currentDensity = db.lookupObject<volVectorField>("totalCurrent");
	volScalarField scalarCurrentDensity = currentDensity.component(0);

	scalar j = (gSum(scalarCurrentDensity)/ncells_)*(gSum(scalarCurrentDensity)/ncells_);

	return j;
}

Foam::scalar Foam::emcModels::current::powerSumMesh() const
{
	const objectRegistry& db = E_.db();
    const volVectorField& tddtE = db.lookupObject<volVectorField>("ddtE");

	volScalarField tpowerSumMesh = meshV_*((mspm_.netChargeFlux() + plasmaConstants::epsilon0*tddtE) & E_);

	return gSum(tpowerSumMesh);
}


void Foam::emcModels::current::correct(dictionary& voltageDict)
{
	if (mode_ == "continuousFrequencyModulated")
	{
		const scalar& tcurrentDensity = currentDensitySum();
		const scalar& tpower = powerSumMesh();

		currentSum_ += (tcurrentDensity+tcurrentDensityOld_)*time_.deltaT().value()*0.5;
		tcurrentDensityOld_ = tcurrentDensity;

		curTimeIndex_ = time_.timeIndex();

		powerSum_ += (tpower+tpowerOld_)*time_.deltaT().value()*0.5;
		tpowerOld_ = tpower;

		if (curTimeIndex_ == 1)
		{
			amplitude_ = initialAmplitude_;
		}

		timeCount_ = timeCount_ + time_.deltaT().value();

		if (timeCount_ >= 1/frequency_)
		{
			scalar rmsCurrent_ = Foam::sqrt(frequency_*currentSum_);
			scalar powerSumAve_ = powerSum_*frequency_;

			powerSum_ = 0.0;
			currentSum_ = 0.0;
			timeCount_ = 0.0;

			scalar amplitudeOld_(amplitude_);

			amplitude_ = amplitudeOld_*(1.0-dampingFactor_*((rmsCurrent_/rms_)-1.0));

			scalar dif = Foam::mag((rmsCurrent_-rms_)/((rmsCurrent_+rms_)/2.0))*100;
			Info << "percentage of difference calculated between desired and" << endl;
			Info << "calculated RMS current is " << dif << endl;

			if (dif<=tolerance_)
			{
				amplitude_=amplitudeOld_;
			}
			else
			{
				if((amplitude_/amplitudeOld_) >= 1.1)
				{
					amplitude_ = 1.1*amplitudeOld_;
				}
				else if((amplitude_/amplitudeOld_) <= 0.9)
				{
					amplitude_ = 0.9*amplitudeOld_;
				}
			}

			if (Pstream::master())
			{
				OFstream& resistanceLogFile = *LogFilePtr_;
				resistanceLogFile << time_.value() << tab;
				resistanceLogFile << amplitude_ << tab;
				resistanceLogFile << rmsCurrent_ << tab;
				resistanceLogFile << powerSumAve_ << endl;
			}
		}

		if (operation_ == "sinusoidal")
		{
			scalar voltageValue = amplitude_*Foam::cos(2.0*M_PI*frequency_*time_.value()) + bias_;
			voltageDict.set("voltage", voltageValue);
		}
		if (operation_ == "cosinusoidal")
		{
			scalar voltageValue = amplitude_*Foam::cos(2.0*M_PI*frequency_*time_.value()) + bias_;
			voltageDict.set("voltage", voltageValue);
		}
		else if (operation_ == "pulsed")
		{
			scalar period = 1/frequency_;
			scalar period_duty = period*dutyCycle_/100;
			const scalar wd_ = w_*Foam::sqrt(1-Foam::sqr(e_));

			scalar n = floor((this->db().time().value())/period);
			scalar pTime_ = db().time().value() - period*n;

			if (waveform_ == "symmetricalBipolar")
			{
				if (this->db().time().value() < period*n+(period/2))
				{
					scalar voltageValue = amplitude_;
					voltageDict.set("voltage", voltageValue);
				}
				else if (this->db().time().value() > period*n+(period/2))
				{
					scalar voltageValue = -amplitude_;
					voltageDict.set("voltage",voltageValue);
				}
			}
			else if (waveform_ == "unipolar")
			{
				if (this->db().time().value() < period*n+period_duty)
				{
					scalar voltageValue = amplitude_;
					voltageDict.set("voltage", voltageValue);
				}
				else if (this->db().time().value() > period*n+period_duty)
				{
					scalar voltageValue = 0;
					voltageDict.set("voltage", voltageValue);
				}
			}
			else if (waveform_ == "bias")
			{
				if (this->db().time().value() < period*n+period_duty)
				{
					scalar voltageValue = amplitude_+bias_;
					voltageDict.set("voltage", voltageValue);
				}
				else if (this->db().time().value() > period*n+period_duty)
				{
					scalar voltageValue = -amplitude_+bias_;
					voltageDict.set("voltage",voltageValue);
				}
			}
			else if (waveform_ == "underDamped")
			{
				if (this->db().time().value() < period*n+period_duty)
				{
					scalar voltageValue = Foam::exp(-e_*w_*pTime_)*(amplitude_*Foam::cos(wd_*pTime_)+(1/wd_+e_*amplitude_-e_*2*Foam::sin(wd_*pTime_))) + bias_;
					voltageDict.set("voltage",voltageValue);
				}
				else if (this->db().time().value() > period*n+period_duty)
				{
					scalar voltageValue = bias_;
					voltageDict.set("voltage",voltageValue);
				}
			}
		}
	}
	else
	{
        FatalErrorIn("emcModels::power::correct(dictionary& voltageDict)")
            << " incorrect mode:  flag 'continuousFrequencyModulated' for"
            << " mode not found"
            << exit(FatalError);
	}
}

bool Foam::emcModels::current::read(const dictionary& electroMagnetics)
{
    emcModel::read(electroMagnetics);
    emcModelCoeffs_.lookup("initialAmplitude") >> initialAmplitude_;
    emcModelCoeffs_.lookup("frequency") >> frequency_;
    emcModelCoeffs_.lookup("bias") >> bias_;
    emcModelCoeffs_.lookup("rms") >> rms_;
    emcModelCoeffs_.lookup("dampingFactor") >> dampingFactor_;
    emcModelCoeffs_.lookup("dutyCycle") >> dutyCycle_;
    emcModelCoeffs_.lookup("naturalFrequency") >> w_;
    emcModelCoeffs_.lookup("dampingRatio") >> e_;
    emcModelCoeffs_.lookup("tolerance") >> tolerance_;

    return true;
}