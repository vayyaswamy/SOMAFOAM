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
				const word mode_;				// flag check
				scalar initialAmplitude_;		// initial amplitude of voltage
				scalar frequency_;				// operational frequency	
				scalar bias_;					// bias in the plasma 
				scalar rms_;					// operational power of the plasma
				scalar dampingFactor_;			// time-advance voltage variation 
				scalar timeCount_;				// using a constant time steps, how many steps needed for cycle of interest
				label curTimeIndex_;			// when period ends, saves ith time index for next iteration

				OFstream* LogFilePtr_;			// stream of data for electrical characterization
				scalar amplitude_;				// amplitude utilized for different waveforms at n time

				scalar rmsCurrent_;				// root-mean-squared current of the system
				scalar ncells_;					// n of cells
				scalar currentSum_;				// summation of the current during the cycle
				const word operation_;			// type of waveform? sinusoidal or pulsed? 
				const word waveform_;			// type of pulsed waveform
				scalar dutyCycle_;				// percentage waveform is on or off
				scalar w_;						// frequency associated with the decay waveform 
				scalar e_;						// damping of the decay waveform
				scalar powerSum_;				// summation of power during a cycle 
				scalar tolerance_;				// percentage of difference parameter between calculated power and user defined power

				scalar tpowerOld_;				// summation of power during a cycle (old value) 
				scalar tcurrentDensityOld_;  	// summation of current during a cycle (old value) 


				volScalarField meshV_;			// mesh volumes

				// Private Member Functions
				scalar currentDensitySum() const;
				scalar powerSumMesh() const;

			    //- Disallow copy construct
				current(const current &);

				//- Disallow default bitwise assignment
				void operator=(const current&);

			public:

			    //- Runtime type information
				TypeName("current");

			    //- Construct from components
				current
				(
					const dictionary& electroMagnetics,
					multiSpeciesPlasmaModel& mspm,
					const volVectorField& E,
					const Time& runTime
				);

			    // Destructor
				virtual ~current();

			    // execution of power controls
				virtual void correct(dictionary& voltageDict);

			    //- Update properties from given dictionary
				virtual bool read(const dictionary& electroMagnetics);
		};
	};
}

#endif 