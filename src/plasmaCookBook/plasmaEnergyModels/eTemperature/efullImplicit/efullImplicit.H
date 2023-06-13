/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::eTemp

Description
    Abstract base class for plasma energy models.

SourceFiles
    eTemp.C

\*---------------------------------------------------------------------------*/

#ifndef efullImplicit_H
#define efullImplicit_H

#include "eTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class efullImplicit Declaration
\*---------------------------------------------------------------------------*/

class efullImplicit
:
    public eTemp
{
    // Private data

        volVectorField eeFlux;

        const word eSpecie;

        const label eIndex_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        efullImplicit(const efullImplicit&);

        //- Disallow default bitwise assignment
        void operator=(const efullImplicit&);


public:

    //- Runtime type information
    TypeName("efullImplicit");

    // Constructors

        //- Construct from components
        efullImplicit
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~efullImplicit()
        {}


    // Member Functions

        virtual Foam::scalar correct(psiChemistryModel& chemistry, const volVectorField& E);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
