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

#ifndef iTemp_H
#define iTemp_H

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
                        Class iTemp Declaration
\*---------------------------------------------------------------------------*/

class iTemp
{
private:

    // Private data

        hsCombustionThermo& thermo_;

        multiSpeciesPlasmaModel& mspm_;

        const volVectorField& E_;

        const fvMesh& mesh_;
    // Private Member Functions

        //- Disallow default bitwise copy construct
        iTemp(const iTemp&);

        //- Disallow default bitwise assignment
        void operator=(const iTemp&);


public:

    //- Runtime type information
    TypeName("iTemp");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            iTemp,
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
        iTemp
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E
        );


    // Selectors

        static autoPtr<iTemp> New
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~iTemp()
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

        virtual void correct(psiChemistryModel& chemistry) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
