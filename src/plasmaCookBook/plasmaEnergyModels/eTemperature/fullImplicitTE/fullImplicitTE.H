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

#ifndef fullImplicitTE_H
#define fullImplicitTE_H

#include "eTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class efullImplicitTE Declaration
\*---------------------------------------------------------------------------*/

class efullImplicitTE
:
    public eTemp
{
    // Private data

        volVectorField eeFlux;

        const word eSpecie;

        const label eIndex_;

        dictionary fieldBounds;

        scalar TeMin;

        scalar TeMax;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        efullImplicitTE(const efullImplicitTE&);

        //- Disallow default bitwise assignment
        void operator=(const efullImplicitTE&);

public:

    //- Runtime type information
    TypeName("efullImplicitTE");

    // Constructors

        //- Construct from components
        efullImplicitTE
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );

    // Destructor

        virtual ~efullImplicitTE()
        {}

    // Member Functions

        virtual void correct(psiChemistryModel& chemistry, const volVectorField& E);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
