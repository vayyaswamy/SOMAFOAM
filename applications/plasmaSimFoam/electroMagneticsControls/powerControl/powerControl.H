/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#ifndef powerControl_H
#define powerControl_H

#include "emcModels.H"

namespace Foam
{
namespace emcModels
{

/*---------------------------------------------------------------------------*\
                          Class power Declaration
\*---------------------------------------------------------------------------*/

class power
:
    public emcModel
{
    // Private data
		const word mode_;
		scalar initialAmplitude_;
		scalar frequency_;
		scalar bias_;
		scalar power_;
		scalar dampingFactor_;
		scalar timeCounter_;
		scalar timeCount_;
		label curTimeIndex_;

		OFstream* powerLogFilePtr_;
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

    // Private Member Functions

		scalar currentDensitySum() const;

        inline scalar powerSumMesh() const;

        //- Disallow copy construct
        power(const power&);

        //- Disallow default bitwise assignment
        void operator=(const power&);


public:

    //- Runtime type information
    TypeName("power");


    // Constructors

        //- Construct from components
        power
        (
			const dictionary& electroMagnetics,
			multiSpeciesPlasmaModel& mspm,
			const volVectorField& E,
			const Time& runTime
        );


    // Destructor

        virtual ~power();

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
