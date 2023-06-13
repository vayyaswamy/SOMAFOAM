/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::iTemp

Description
    Abstract base class for plasma energy models.

SourceFiles
    iTemp.C

\*---------------------------------------------------------------------------*/

#ifndef ifullImplicit_H
#define ifullImplicit_H

#include "iTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ifullImplicit Declaration
\*---------------------------------------------------------------------------*/

class ifullImplicit
:
    public iTemp
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ifullImplicit(const ifullImplicit&);

        //- Disallow default bitwise assignment
        void operator=(const ifullImplicit&);


public:

    //- Runtime type information
    TypeName("ifullImplicit");

    // Constructors

        //- Construct from components
        ifullImplicit
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~ifullImplicit()
        {}


    // Member Functions

        virtual void correct(psiChemistryModel& chemistry);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
