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
    dampingFactor_(readScalar(emcModelCoeffs_.lookup("dampingFactor"))),
    ncells_(0.0),
    rmsCount_(1.0),
	currentSum_(0.0),
	operation_(emcModelCoeffs_.lookup("waveform")),
	waveform_(emcModelCoeffs_.lookup("operation")),
	dutyCycle_(readScalar(emcModelCoeffs_.lookup("dutyCycle"))),
	w_(readScalar(emcModelCoeffs_.lookup("naturalFrequency"))),
	e_(readScalar(emcModelCoeffs_.lookup("dampingRatio"))),
	powerSum_(0.0),
	tolerance_(readScalar(emcModelCoeffs_.lookup("tolerance"))),
	timeCounter_(0.0),
	timeCount_(0.0),
	curTimeIndex_(time_.timeIndex()),
	powerLogFilePtr_(NULL),
    amplitude_(0.0),
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
        powerLogFilePtr_ = new OFstream(fileName("power_voltage_logfile"));
        OFstream& powerLogFile = *powerLogFilePtr_;
    	powerLogFile << "time" << tab;
    	powerLogFile << "input" << tab;
    	powerLogFile << "RMScurrent" << tab;
    	powerLogFile << "power" << endl;
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

Foam::scalar Foam::emcModels::power::currentDensitySum() const
{
	const objectRegistry& db = E_.db();
	const volVectorField& currentDensity = db.lookupObject<volVectorField>("totalCurrent");
	volScalarField scalarCurrentDensity = currentDensity.component(0);

	scalar j = (gSum(scalarCurrentDensity)/ncells_)*(gSum(scalarCurrentDensity)/ncells_);

	return j;
}

void Foam::emcModels::power::correct(dictionary& voltageDict)
{
	if (mode_ == "continuousFrequencyModulated")
	{
		const scalar& tpower = powerSumMesh();
		const scalar& tcurrentDensity = currentDensitySum();

		powerSum_ += tpower;
		currentSum_ += tcurrentDensity;

		curTimeIndex_ = time_.timeIndex();

		if(curTimeIndex_ == 1)
		{
			amplitude_ = initialAmplitude_;
		}

		timeCount_ = 1/frequency_/time_.deltaT().value();

		if(curTimeIndex_ - timeCounter_ >= timeCount_)
		{
			timeCounter_ =  curTimeIndex_;

			scalar rmsCurrent_ = Foam::sqrt(currentSum_/rmsCount_);
			scalar powerSumAve_ = powerSum_/timeCount_;

			powerSum_ = 0.0;
			currentSum_ = 0.0;
			rmsCount_ = 0.0;

			scalar amplitudeOld_(amplitude_);

			amplitude_ = amplitudeOld_*(1.0-dampingFactor_*((powerSumAve_/power_)-1.0));

			scalar dif = Foam::mag((powerSumAve_-power_)/((powerSumAve_+power_)/2.0))*100;
			Info << "percentage of difference calculated between desired and" << endl;
			Info << "calculated average power is " << dif << endl;

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
   				OFstream& resistanceLogFile = *powerLogFilePtr_;
				resistanceLogFile << time_.value() << tab;
				resistanceLogFile << amplitude_ << tab;
				resistanceLogFile << rmsCurrent_ << tab;
				resistanceLogFile << powerSumAve_ << endl;
			}
		}
		else
		{
			rmsCount_ = rmsCount_ + 1;
		}

		if (operation_ == "sinusoidal")
		{
			scalar voltageValue = amplitude_*Foam::sin(2.0*M_PI*frequency_*time_.value()) + bias_;
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
    emcModelCoeffs_.lookup("dampingFactor") >> dampingFactor_;
    emcModelCoeffs_.lookup("dutyCycle") >> dutyCycle_;
    emcModelCoeffs_.lookup("naturalFrequency") >> w_;
    emcModelCoeffs_.lookup("dampingRatio") >> e_;
    emcModelCoeffs_.lookup("tolerance") >> tolerance_;

    return true;
}


// ************************************************************************* //
