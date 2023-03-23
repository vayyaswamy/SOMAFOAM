#ifndef one_H
#define one_H

#include "spatioTemporal.H"

namespace Foam
{
	namespace spatioTemporals
	{
		class one
		:
		public spatioTemporal
		{
			private:

				const scalar frequency_;					// frequency of operation
				const scalar nCycles_;						// number of cycles to record

				scalar startWriting_;						// time recording of data starts

				PtrList<volScalarField> Nspecie_;			// list recording species number density

				const word temporal_[7];					// array containing data sets types IO parameters
				const word units_[7];						// array containing data sets types units
				const int streamSize_;						// size of the IO stream 
				List<autoPtr<Foam::OFstream> > streamers_;	// data streams to IO into output files	
				const vectorField& xCoord_;					// x-coordinates centroids 

				const volScalarField& Phi_;					// potential across gap
				const volScalarField& Te_;					// electron energy across gap
				const volVectorField& ddtE_;				// electric field 
				const volVectorField& currentDensity_;		// current density 
				const volScalarField& eDivergence_;			// electron's sources and sinks

				// disallow copy construct 
				one(const one&);

				// disallow default bitwise assigment 
				void operator=(const one&);

			public:

				TypeName("oneD");

				// constructor
				one
				(
					const dictionary& spatioTemporals,
					const dictionary& electroMagnetics,
					multiSpeciesPlasmaModel& mspm,
					const volVectorField& E,
					const Time& runTime,
					const fvMesh& mesh
				);

				// destroyer
				virtual ~one();

				// write functions processing a volScalarField
				virtual void write(const volScalarField& input, Ostream& logFile);

				// write functions processing a volVectorField
				virtual void write(const volVectorField& input, Ostream& logFile);

				// execution of IO algorithm
				virtual void correct();	
		};
	}
}

#endif
