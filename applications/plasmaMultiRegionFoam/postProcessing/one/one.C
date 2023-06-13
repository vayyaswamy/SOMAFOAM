#include "one.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
	namespace spatioTemporals
	{
		defineTypeNameAndDebug(one, 0);
		addToRunTimeSelectionTable(spatioTemporal, one, dictionary);
	};
};

Foam::spatioTemporals::one::one
(
	const dictionary& spatioTemporals,
	const dictionary& electroMagnetics,
	multiSpeciesPlasmaModel& mspm,
	const volVectorField& E,
	const Time& runTime,
	const fvMesh& mesh
)
:
	spatioTemporal(spatioTemporals,electroMagnetics,mspm,E,runTime,mesh),
	frequency_(readScalar(spatioTemporalcoeffs_.lookup("frequency"))),
	nCycles_(readScalar(spatioTemporalcoeffs_.lookup("nCycles"))),
	startWriting_(time_.endTime().value()-(1/frequency_)*nCycles_),
	Nspecie_(mspm_.species().size()),
	temporal_({"Phi","Te","Ex","power","currentDensity","conductionCurrent","displacementCurrent"}),
	units_({"V", "K", "V/m", "W/m3", "J", "J", "J"}),
	streamSize_(sizeof(temporal_)/sizeof(temporal_[0])+mspm_.species().size()+1),
	streamers_(streamSize_),
	xCoord_(mesh_.C().internalField()),
	Phi_(mesh_.lookupObject<volScalarField>("Phi")),
	Te_(mesh_.lookupObject<volScalarField>("Te")),
	ddtE_(mesh_.lookupObject<volVectorField>("ddtE")),
	currentDensity_(mesh_.lookupObject<volVectorField>("totalCurrent")),
	eDivergence_(mesh_.lookupObject<volScalarField>("electronDiv"))
	{
		// checks if number of cycles can be recorded with the given end time of simulation
		if (startWriting_ <= 0)
		{
			FatalErrorIn
			(
				"Foam::spatioTemporals::one::one"
			) << "Number of cycles " << nCycles_ << " for a operating frequency of " << frequency_ << " requires more time than the time allocated for the simulation" << endl << exit(FatalError);
		}

		// checks if mesh is a 1D mesh
		if (mesh_.nSolutionD() != 1)
		{
			FatalErrorIn
			(
				"Foam::spatioTemporals::one::one"
			) << "dimensions of the mesh are not 1D" << endl << exit(FatalError);
		}

		// set all the folders for species in the system, and the IO files
		forAll(mspm_.species(),i)
		{
			mkDir(time_.path()/"spatioTemporal"/mspm_.species()[i]);

			streamers_[i].set
			(
				new Foam::OFstream
				(
					fileName(time_.path()/"spatioTemporal"/mspm_.species()[i]/mspm_.species()[i]+".dat")
				)
			);

			Ostream& logFile = streamers_[i]()();
			logFile << "time" << tab << "position" << tab << "density" << endl;
		}

		// set up folders for the other parameters in the system 
		for (int i=0; i<(sizeof(temporal_)/sizeof(temporal_[0])); i++)
		{
			mkDir(time_.path()/"spatioTemporal"/temporal_[i]);
		}

		mkDir(time_.path()/"spatioTemporal"/"electrons-1");

		// set all the folders for other parameters probed in the system, and the IO files
		int ii=0;
		for (int i=mspm_.species().size(); i<streamSize_; i++)
		{
			if (i == streamSize_-1)
			{
				streamers_[i].set
				(
					new Foam::OFstream
					(
						fileName(time_.path()/"spatioTemporal"/"electrons-1"/"electrons-1"+".dat")
					)
				);

				Ostream& logFile = streamers_[i]()();
				logFile << "time" << tab << "position" << tab << "1/m3s-1" << endl;
			}
			else
			{
				streamers_[i].set
				(
					new Foam::OFstream
					(
						fileName(time_.path()/"spatioTemporal"/temporal_[ii]/temporal_[ii]+".dat")
					)
				);

				Ostream& logFile = streamers_[i]()();
				logFile << "time" << tab << "position" << tab << units_[ii] << endl;
				ii = ii + 1;
			}

		}
	}

Foam::spatioTemporals::one::~one()
{

}

void Foam::spatioTemporals::one::write(const volScalarField& input, Ostream& logFile)
{
	forAll(input, celli)
	{
		logFile << time_.value() << tab
		<< xCoord_[celli].component(0) << tab
		<< input[celli] << endl;
	}
}

void Foam::spatioTemporals::one::write(const volVectorField& input, Ostream& logFile)
{
	forAll(input, celli)
	{
		logFile << time_.value() << tab
		<< xCoord_[celli].component(0) << tab
		<< input[celli] << endl;
	}
}

void Foam::spatioTemporals::one::correct()
{
	if (time_.value() >= startWriting_)
	{

		Info << "spatio-temporal IO algorithm in process" << endl;

		volVectorField netChargeFlux_
		(
			IOobject
			(
				"netChargeFlux_",
				mesh_.time().timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mspm_.netChargeFlux()
		);

		volScalarField power_
		(
			IOobject
			(
				"power_",
				time_.timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			((netChargeFlux_.component(0)+plasmaConstants::epsilon0*ddtE_.component(0))*E_.component(0))
		);

		volScalarField disCurrent_
		(
			IOobject
			(
				"disCurrent_",
				time_.timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			plasmaConstants::epsilon0*ddtE_.component(0)
		);

		volScalarField conCurrent_
		(
			IOobject
			(
				"conCurrent_",
				time_.timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			netChargeFlux_.component(0)
		);

		volScalarField divElectron 
		(
			IOobject
			(
				"divElectron",
				time_.timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			eDivergence_
		);

		// write species number density into corresponding IO files
		forAll(mspm_.species(),i)
		{
			Nspecie_.set
			(
				i,
				volScalarField
				(
					IOobject
					(
						mspm_.species()[i]+"_",
						time_.timeName(),
						mesh_,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mspm_.N(i)
				)
			);
			write(Nspecie_[i],streamers_[i]()());
		}

		// index for data stream after writing species IO
		int j = mspm_.species().size();

		// write measured parameters into corresponding IO files
		write(Phi_, streamers_[j]()());
		write(Te_, streamers_[j+1]()());
		write(E_,streamers_[j+2]()());
		write(power_,streamers_[j+3]()());
		write(currentDensity_, streamers_[j+4]()());
		write(conCurrent_, streamers_[j+5]()());
		write(disCurrent_, streamers_[j+6]()());
		write(divElectron, streamers_[j+7]()());
	}
}
