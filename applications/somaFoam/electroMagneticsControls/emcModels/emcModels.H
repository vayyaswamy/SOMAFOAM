/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::emcModel

Description
    Base-class for electromagnetic controls.

SourceFiles
    emcModel.C

\*---------------------------------------------------------------------------*/

#ifndef emcModels_H
#define emcModels_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "multiSpeciesPlasmaModel.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class emcModel Declaration
\*---------------------------------------------------------------------------*/

class emcModel
:
    public regIOobject
{

protected:

    // Protected data

        dictionary emcModelCoeffs_;

        multiSpeciesPlasmaModel& mspm_;
        const volVectorField& E_;
		const Time& time_;

private:

    // Private Member Functions

        //- Disallow copy construct
        emcModel(const emcModel&);

        //- Disallow default bitwise assignment
        void operator=(const emcModel&);


public:

    //- Runtime type information
    TypeName("emcModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            emcModel,
            dictionary,
            (
				const dictionary& electroMagnetics,
				multiSpeciesPlasmaModel& mspm,
				const volVectorField& E,
				const Time& runTime
            ),
            (
                electroMagnetics,
                mspm,
                E,
				runTime
            )
        );


    // Selectors

        static autoPtr<emcModel> New
        (
			const dictionary& electroMagnetics,
			multiSpeciesPlasmaModel& mspm,
			const volVectorField& E,
			const Time& runTime
        );


    // Constructors

        //- Construct from components
        emcModel
        (
			const dictionary& electroMagnetics,
			multiSpeciesPlasmaModel& mspm,
			const volVectorField& E,
			const Time& runTime
        );


    // Destructor

        virtual ~emcModel();


    // Member Functions

        virtual void correct(dictionary& voltageDict) = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& PDRProperties) = 0;

        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
