#ifndef spatioTemporal_H
#define spatioTemporal_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "multiSpeciesPlasmaModel.H"
#include <math.h>

namespace Foam
{
	class spatioTemporal
	:
	public regIOobject
	{
		protected:

			dictionary spatioTemporalcoeffs_;
			dictionary emcModelCoeffs_;
			multiSpeciesPlasmaModel& mspm_;
			const volVectorField& E_;
			const Time& time_;
	        const fvMesh& mesh_;

		private:

			//- Disallow copy construct
			spatioTemporal(const spatioTemporal&);

			// disallow default bitwise assigment
			void operator=(const spatioTemporal&);

		public:

			// Runtime type information
			TypeName("spatioTemporal");

			declareRunTimeSelectionTable
			(
				autoPtr,
				spatioTemporal,
				dictionary,
				(
					const dictionary& spatioTemporals,
					const dictionary& electroMagnetics,
					multiSpeciesPlasmaModel& mspm,
					const volVectorField& E,
					const Time& runTime,
					const fvMesh& mesh
				),
				(
					spatioTemporals,
					electroMagnetics,
					mspm,
					E,
					runTime,
					mesh
				)
			);

			// selector
			static autoPtr<spatioTemporal> New
			(
				const dictionary& spatioTemporals,
				const dictionary& electroMagnetics,
				multiSpeciesPlasmaModel& mspm,
				const volVectorField& E,
				const Time& runTime,
				const fvMesh& mesh
			);

			// constructor
			spatioTemporal
			(
				const dictionary& spatioTemporals,
				const dictionary& electroMagnetics,
				multiSpeciesPlasmaModel& mspm,
				const volVectorField& E,
				const Time& runTime,
				const fvMesh& mesh
			);

			// destructor
			virtual ~spatioTemporal();

			virtual void correct() = 0;

	        virtual bool writeData(Ostream&) const
	        {
	            return true;
	        }
	};
}

#endif