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


    // Protected member functions

        //- Drive value at time t
        scalar driveValue(const scalar t) const;

        //- Set Te from the drive at the current time
        void setDrivenTe();

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
