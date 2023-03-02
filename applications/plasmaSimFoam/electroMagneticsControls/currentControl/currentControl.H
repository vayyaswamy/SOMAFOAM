#ifndef currentControl_H
#define currentControl_H

#include "emcModels.H"

namespace Foam
{
	namespace emcModels
	{
		class current
		: 
		public emcModel
		{
			private:
				const word mode_;
				scalar initialAmplitude_;
				scalar frequency_;
				scalar bias_;
				scalar rms_;
				scalar dampingFactor_;
				scalar timeCounter_;
				scalar timeCount_;
				label curTimeIndex_;

				OFstream* LogFilePtr_;
				scalar amplitude_;
				scalar rmsCurrent_;
				scalar ncells_;
				scalar rmsCount_;
				scalar currentSum_;
				const word operation_;
				const word waveform_;
				scalar dutyCycle_;
				scalar w_;
				scalar e_;
				scalar powerSum_;
				scalar tolerance_;
				volScalarField meshV_;


				scalar currentDensitySum() const;

				scalar powerSumMesh() const;

				current(const current &);

				void operator=(const current&);

			public:

				TypeName("current");

				current
				(
					const dictionary& electroMagnetics,
					multiSpeciesPlasmaModel& mspm,
					const volVectorField& E,
					const Time& runTime
				);

				virtual ~current();

				virtual void correct(dictionary& voltageDict);

				virtual bool read(const dictionary& electroMagnetics);
		};
	};
}

#endif 