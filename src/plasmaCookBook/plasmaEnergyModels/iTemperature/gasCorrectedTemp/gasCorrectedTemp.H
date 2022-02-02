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

#ifndef gasCorrectedTemp_H
#define gasCorrectedTemp_H

#include "iTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gasCorrectedTemp Declaration
\*---------------------------------------------------------------------------*/

class gasCorrectedTemp
:
    public iTemp
{
    // Private data

        const word ionSpecie;

        const label iIndex_;

        const word backgroundGas;

        const label bIndex_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gasCorrectedTemp(const gasCorrectedTemp&);

        //- Disallow default bitwise assignment
        void operator=(const gasCorrectedTemp&);


public:

    //- Runtime type information
    TypeName("gasCorrectedTemp");

    // Constructors

        //- Construct from components
        gasCorrectedTemp
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~gasCorrectedTemp()
        {}


    // Member Functions

        virtual void correct(psiChemistryModel& chemistry);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
