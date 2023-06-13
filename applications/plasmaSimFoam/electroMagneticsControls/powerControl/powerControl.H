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
	const word mode_;				// flag check
	scalar initialAmplitude_;		// initial amplitude of voltage
	scalar frequency_;				// operational frequency
	scalar bias_;					// bias in the plasma 
	scalar power_;					// operational power of the plasma
	scalar dampingFactor_;			// time-advance voltage variation 
	scalar timeCount_;				// counts dt in system
	label curTimeIndex_;			// initialize initial voltage set by the boundary condition

	OFstream* powerLogFilePtr_;		// stream of data for electrical characterization
	scalar amplitude_;				// amplitude utilized for different waveforms at n time

	scalar ncells_;					// n of cells
	scalar currentSum_;				// summation of the current during the cycle
	const word operation_;			// type of waveform? sinusoidal or pulsed? 
	const word waveform_;			// type of pulsed waveform
	scalar dutyCycle_;				// percentage waveform is on or off
	scalar w_;						// frequency associated with the decay waveform 
	scalar e_;						// damping of the decay waveform
	scalar powerSum_;				// summation of power during a cycle 
	scalar tpowerOld_;				// summation of power during a cycle (old value) 
	scalar tcurrentDensityOld_;  	// summation of current during a cycle (old value) 

	scalar tolerance_;				// percentage of difference parameter between calculated power and user defined power

    volScalarField meshV_;			// mesh volumes

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

    // execution of power controls
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
