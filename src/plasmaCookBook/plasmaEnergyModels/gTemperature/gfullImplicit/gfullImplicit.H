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

#ifndef gfullImplicit_H
#define gfullImplicit_H

#include "gTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gfullImplicit Declaration
\*---------------------------------------------------------------------------*/

class gfullImplicit
:
    public gTemp
{
    // Private data

        const word tavesource;

        const scalar cycleavevalue;

        volScalarField gasTempSource;

        volScalarField gasTempSourcetemp;

        scalar timeCount;
        
        scalar timeCountn;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gfullImplicit(const gfullImplicit&);

        //- Disallow default bitwise assignment
        void operator=(const gfullImplicit&);


public:

    //- Runtime type information
    TypeName("gfullImplicit");

    // Constructors

        //- Construct from components
        gfullImplicit
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~gfullImplicit()
        {}


    // Member Functions

        virtual void correct(psiChemistryModel& chemistry, const volVectorField& E);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
