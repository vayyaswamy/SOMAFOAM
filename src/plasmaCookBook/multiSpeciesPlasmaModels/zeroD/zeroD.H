/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Class
    Foam::zeroD

Description
    Zero-dimensional (global) model: the species change by the chemistry
    alone, dN/dt = RR, on a mesh of one cell whose faces are all empty. Use
    it with plasmaChemistryModel temporal, so that RR is the rate averaged
    over the time step by the ODE solver of the chemistry, and with
    poissonEquationSolver none in constant/electroMagnetics.

    The discharge is driven by the optional dictionary zeroDDrive in
    constant/plasmaProperties

        zeroDDrive
        {
            mode    reducedField;   // reducedField, electronTemperature or none
            table   ((0 100) (1e-8 100) (1.1e-8 0));   // (time value)
        }

    or with "value 100;" in place of the table. The table is interpolated
    linearly in time and held at its end values outside it.

    reducedField: E/N in Td. The field E = (E/N) N along x sets the electron
    flux mu_e N_e z_e E, whose work heats the electrons in the electron
    energy equation. The electron mobility must be eTemp or constant.

    electronTemperature: Te in K, set at the start of every step; the
    electron energy equation is not solved.

    power: power absorbed by the electrons in W. The field of each step is
    set from the electron density and mobility at its start so that the
    electrons absorb this power in the volume of zeroDWall (or "volume"
    in zeroDDrive), e Ne mu_e E^2 V = P; the electron flux is then set so
    that the work J.E of the step is exactly P/V. The field is limited to
    maxReducedField (Td, 3000 by default) when there are few electrons.
    The electron mobility must be eTemp or constant.

    The optional dictionary zeroDWall adds the losses to the walls of a
    cylinder of radius R and length L, as in global models (Thorsteinsson
    and Gudmundsson, Plasma Sources Sci. Technol. 19, 015001 (2010),
    equations (12)-(19), after Kim et al., J. Vac. Sci. Technol. A 24,
    2025 (2006)):

        zeroDWall
        {
            radius      0.0508;     // [m]
            length      0.0254;     // [m]
            ionNeutralCrossSection 1e-18;   // [m2], lambda_i = 1/(N sigma)
            ionIonRecombination    4e-13;   // [m3/s], k_rec of n_*
            ions        ((CF3p1 CF3) (O2p1 O2));    // (ion neutral)
            neutrals
            {
                F   { gamma 0.02;  DN 1e21; pool yes; products ((F2 0.5)); }
                CF3 { gamma 0.017; DN 1e21; fromPool 1; products ((CF4 1)); }
                O   { gamma 0.4;   DN 1e21; products ((O2 0.5)); }
            }
        }

    Positive ions are lost at u_B0 A_eff/V with u_B0 = (kB Te/m_i)^1/2 and
    A_eff = 2 pi (R^2 h_L + R L h_R); the edge-to-centre ratios h_L and h_R
    take the electronegativity alpha = n-/ne, with alpha_0 = 1.5 alpha and
    the positive and negative ions at the gas temperature. Each ion returns
    as its neutral, and one electron is lost with it; negative ions are not
    lost. The electrons lost carry 2 kB Te each out of the electron energy
    (the field electronWallEnergyLoss, W/m3, read by the energy equation).
    A neutral sticks at the rate
        k = [Lambda^2/D + 2 V (2 - gamma)/(A v gamma)]^-1,
    Lambda^-2 = (pi/L)^2 + (2.405/R)^2, D = DN/N, v = (8 kB Tg/(pi m))^1/2,
    and returns as its products. The atoms of the species marked "pool"
    are shared with the species that need "fromPool" of them to form their
    products (CF3 + F -> CF4 at the wall); when the pool is short, that
    fraction of each such species is reflected unchanged, and what is left
    of the pool returns as the products of the pool species. All wall
    losses of a step are exact for the densities of its start,
    N (1 - exp(-k dt)).

SourceFiles
    zeroD.C

\*---------------------------------------------------------------------------*/

#ifndef zeroD_H
#define zeroD_H
#include "fvm.H"
#include "Tuple2.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class zeroD Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class zeroD
:
    public multiSpeciesPlasmaModel
{

protected:

    // Protected data

        const PtrList<ThermoType>& speciesThermo_;

        //- Drive: none, reducedField or electronTemperature
        word driveMode_;

        //- Drive value versus time
        List<Tuple2<scalar, scalar> > driveTable_;

        //- Volume [m3] for the power drive
        scalar volume_;

        //- Largest reduced field [Td] of the power drive
        scalar maxReducedField_;

        //- Field of the power drive [V/m], set by setDrivenField
        mutable scalarField powerField_;

        // Wall losses (zeroDWall)

            bool wall_;

            scalar radius_;

            scalar length_;

            scalar sigmaIn_;

            scalar kRec_;

            //- Positive ions lost at the wall and their neutrals
            labelList wallIon_;

            labelList wallIonProduct_;

            //- Neutrals lost at the wall
            labelList wallNeutral_;

            scalarList wallGamma_;

            scalarList wallDN_;

            //- Pool atoms needed per event (0: none)
            scalarList wallFromPool_;

            //- Products (species, number per lost particle)
            List<List<Tuple2<label, scalar> > > wallProducts_;

            //- Position of the pool species in wallNeutral_ (-1: none)
            label wallPool_;

            //- Electron energy carried to the walls [W/m3]
            autoPtr<volScalarField> electronWallEnergyLoss_;


    // Protected member functions

        //- Drive value at time t
        scalar driveValue(const scalar t) const;

        //- Set Te from the drive at the current time
        void setDrivenTe();

        //- Read zeroDWall
        void readWall(const dictionary& dict);

        //- Losses to the walls over a step of length dt, applied to the
        //  densities and mass fractions of the species
        void wallLosses(const scalar dt);

public:

        //- Runtime type information
        TypeName("zeroD");

    // Constructors

        //- Construct from components
        zeroD
        (
            hsCombustionThermo& thermo
        );

        //- Destructor
        ~zeroD()
        {}

    // Member Functions

        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }

        virtual inline scalar correct
        (
            psiChemistryModel& chemistry,
            const volVectorField& E,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );

        //- Read the species input, then check and apply the drive
        virtual void input();

        virtual void setDrivenField(volVectorField& E) const;

        virtual bool fixesTe() const
        {
            return driveMode_ == "electronTemperature";
        }

        //- Read dictionary
        virtual bool read();
  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "zeroD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
