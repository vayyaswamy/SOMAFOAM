/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::eEnergyDensity

Description
    Electron energy equation for the energy density n_eps = 3/2 n_e kB Te,

        d(n_eps)/dt + div(Gamma_eps) = -e Gamma_e.E - S

    with the energy flux in drift-diffusion form (Hagelaar and Pitchford,
    Plasma Sources Sci. Technol. 14 (2005) 722)

        Gamma_eps = 5/3 (F_e n_eps - D_e grad(n_eps))

    F_e and D_e being the drift velocity and diffusivity of the electrons
    (the heat conductivity is then 5/2 kB n_e D_e). The energy flux is
    discretised as the electron flux (fvSchemes, or Scharfetter-Gummel), and
    the energy is conserved by the discretisation. The energy loss S of the
    reactions is linearised in n_eps through Te = n_eps/(3/2 kB n_e); the
    work of the field is that of the electron face flux. The wall value of
    n_eps is 3/2 kB n_e Te with the wall temperature of the Te boundary
    condition, so the wall energy model is that of the Te condition. With a
    gas flow, the energy is carried by the gas.

    Selected in constant/plasmaProperties:
    \verbatim
    energyModel
    {
        electron eEnergyDensity;
        ...
    }
    \endverbatim

    Passes per time step and their tolerance: electronTemperature and
    tolerance of innerIterations, as efullImplicit.

SourceFiles
    eEnergyDensity.C

\*---------------------------------------------------------------------------*/

#ifndef eEnergyDensity_H
#define eEnergyDensity_H

#include "eTemp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eEnergyDensity Declaration
\*---------------------------------------------------------------------------*/

class eEnergyDensity
:
    public eTemp
{
    // Private data

        //- Energy flux 5/2 kB Te J_e (for the time step limit)
        volVectorField eeFlux;

        const word eSpecie;

        const label eIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        eEnergyDensity(const eEnergyDensity&);

        //- Disallow default bitwise assignment
        void operator=(const eEnergyDensity&);


public:

    //- Runtime type information
    TypeName("eEnergyDensity");

    // Constructors

        //- Construct from components
        eEnergyDensity
        (
            hsCombustionThermo& thermo,
            multiSpeciesPlasmaModel& mspm,
            const volVectorField& E,
            const dictionary& dict
        );


    // Destructor

        virtual ~eEnergyDensity()
        {}


    // Member Functions

        virtual Foam::scalar correct
        (
            psiChemistryModel& chemistry,
            const volVectorField& E
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
