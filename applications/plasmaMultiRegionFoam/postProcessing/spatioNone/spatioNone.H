#ifndef spatioNone_H
#define spatioNone_H

#include "spatioTemporal.H"

namespace Foam
{
	namespace spatioTemporals
	{
		class spatioNone
		:
		public spatioTemporal
		{
			private:

				// disallow copy construct 
				spatioNone(const spatioNone&);

				// disallow default bitwise assigment
				void operator=(const spatioNone&);

			public:

				TypeName("noneD");

				spatioNone
				(
					const dictionary& spatioTemporals,
					const dictionary& electroMagnetics,
					multiSpeciesPlasmaModel& mspm,
					const volVectorField& E,
					const Time& runTime,
					const fvMesh& mesh
				);

				virtual ~spatioNone();

				virtual void correct();
		};
	}
}

#endif