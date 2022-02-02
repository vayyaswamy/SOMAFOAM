/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::gTemp

Description
    Abstract base class for plasma energy models.

SourceFiles
    gTemp.C

\*---------------------------------------------------------------------------*/

#ifndef none_H
#define none_H

#include "gTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public gTemp
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
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
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~none()
        {}


    // Member Functions

        virtual void correct(psiChemistryModel& chemistry, const volVectorField& E);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
