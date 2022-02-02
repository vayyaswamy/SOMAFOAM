/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#ifndef none_H
#define none_H

#include "emcModels.H"

namespace Foam
{
namespace emcModels
{

/*---------------------------------------------------------------------------*\
                          Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public emcModel
{
    // Private data

    // Private Member Functions

        //- Disallow copy construct
        none(const none&);

        //- Disallow default bitwise assignment
        void operator=(const none&);


public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        none
        (
			const dictionary& electroMagnetics,
			multiSpeciesPlasmaModel& mspm,
			const volVectorField& E,
			const Time& runTime
        );


    // Destructor

        virtual ~none();


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
