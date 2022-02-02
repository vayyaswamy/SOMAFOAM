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

#ifndef eTemp_H
#define eTemp_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "multiSpeciesPlasmaModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class eTemp Declaration
\*---------------------------------------------------------------------------*/

class eTemp
{
private:

    // Private data

        hsCombustionThermo& thermo_;

        multiSpeciesPlasmaModel& mspm_;

        const volVectorField& E_;

        const fvMesh& mesh_;

		Time& runTime_;

		Switch restartcapable;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        eTemp(const eTemp&);

        //- Disallow default bitwise assignment
        void operator=(const eTemp&);

public:

    //- Runtime type information
    TypeName("eTemp");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            eTemp,
            dictionary,
            (
                hsCombustionThermo& thermo,
                multiSpeciesPlasmaModel& mspm,
                const volVectorField& E,
                const dictionary& dict
            ),
            (thermo, mspm, E, dict)
        );

    // Constructors

        //- Construct from components
        eTemp
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E
        );

    // Selectors

        static autoPtr<eTemp> New
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );

    // Destructor

        virtual ~eTemp()
        {}

    // Member Functions

        hsCombustionThermo& thermo()
        {
            return thermo_;
        }

        multiSpeciesPlasmaModel& mspm()
        {
            return mspm_;
        }

        const volVectorField& E() const
        {
            return E_;
        }

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        Time& runTime()
        {
            return runTime_;
        }

        Switch restartCapable()
        {
            return restartcapable;
        }

        virtual void correct(psiChemistryModel& chemistry, const volVectorField& E) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
