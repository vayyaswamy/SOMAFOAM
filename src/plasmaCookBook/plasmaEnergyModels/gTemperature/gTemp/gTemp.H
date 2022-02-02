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

#ifndef gTemp_H
#define gTemp_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "multiSpeciesPlasmaModel.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class gTemp Declaration
\*---------------------------------------------------------------------------*/

class gTemp
{
private:

    // Private data

        hsCombustionThermo& thermo_;

        multiSpeciesPlasmaModel& mspm_;

        const volVectorField& E_;

        const fvMesh& mesh_;

		Time& runTime_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gTemp(const gTemp&);

        //- Disallow default bitwise assignment
        void operator=(const gTemp&);


public:

    //- Runtime type information
    TypeName("gTemp");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            gTemp,
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
        gTemp
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E
        );


    // Selectors

        static autoPtr<gTemp> New
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~gTemp()
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

        virtual void correct(psiChemistryModel& chemistry, const volVectorField& E) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
