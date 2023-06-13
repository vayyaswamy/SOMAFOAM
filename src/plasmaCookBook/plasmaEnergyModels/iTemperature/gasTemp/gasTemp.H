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

#ifndef gasTemp_H
#define gasTemp_H

#include "iTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gasTemp Declaration
\*---------------------------------------------------------------------------*/

class gasTemp
:
    public iTemp
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gasTemp(const gasTemp&);

        //- Disallow default bitwise assignment
        void operator=(const gasTemp&);


public:

    //- Runtime type information
    TypeName("gasTemp");

    // Constructors

        //- Construct from components
        gasTemp
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~gasTemp()
        {}


    // Member Functions

        virtual void correct(psiChemistryModel& chemistry);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
