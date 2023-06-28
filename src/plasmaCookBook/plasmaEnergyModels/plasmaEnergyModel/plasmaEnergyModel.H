/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::plasmaEnergyModel

Description
    Abstract base class for plasma energy models.

SourceFiles
    plasmaEnergyModel.C

\*---------------------------------------------------------------------------*/


#ifndef plasmaEnergyModel_H
#define plasmaEnergyModel_H

#include "multiSpeciesPlasmaModel.H"

#include "eTemp.H"
#include "iTemp.H"
#include "gTemp.H"

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plasmaEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaEnergyModel
:
    public IOdictionary
{
private:

    // Private data
        autoPtr<eTemp> eTempPtr_;
        autoPtr<iTemp> iTempPtr_;
        autoPtr<gTemp> gTempPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        plasmaEnergyModel(const plasmaEnergyModel&);

        //- Disallow default bitwise assignment
        void operator=(const plasmaEnergyModel&);


public:

    //- Runtime type information
    TypeName("plasmaEnergyModel");


    // Constructors

        //- Construct from components
        plasmaEnergyModel
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E
        );


    // Destructor

        virtual ~plasmaEnergyModel()
        {}


    // Member Functions

        virtual Foam::scalar ecorrect(psiChemistryModel& chemistry, const volVectorField& E);
        virtual void icorrect(psiChemistryModel& chemistry);
        virtual void gcorrect(psiChemistryModel& chemistry, const volVectorField& E);

        //- Read plasmaProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
