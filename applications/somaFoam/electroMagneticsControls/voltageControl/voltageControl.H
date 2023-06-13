/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#ifndef voltageControl_H
#define voltageControl_H

#include "emcModels.H"

namespace Foam
{
namespace emcModels
{

/*---------------------------------------------------------------------------*\
                          Class voltage Declaration
\*---------------------------------------------------------------------------*/

class voltage
:
    public emcModel
{
    // Private data

		const word mode_;		
		dimensionedScalar amplitude_;
		dimensionedScalar frequency_;
		dimensionedScalar bias_;

    // Private Member Functions

        //- Disallow copy construct
        voltage(const voltage&);

        //- Disallow default bitwise assignment
        void operator=(const voltage&);


public:

    //- Runtime type information
    TypeName("voltage");

    // Constructors

        //- Construct from components
        voltage
        (
			const dictionary& electroMagnetics,
			multiSpeciesPlasmaModel& mspm,
			const volVectorField& E,
			const Time& runTime
        );


    // Destructor

        virtual ~voltage();


    // Member Functions

        virtual void correct(dictionary& voltageDict);

        //- Update properties from given dictionary
        virtual bool read(const dictionary& electroMagnetics);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
