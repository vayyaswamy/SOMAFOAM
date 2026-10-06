/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::rfPowerElectrode

Description
    Potential of an RF electrode driven at a given mean power through a
    blocking capacitor:

        V(t) = A sin(2 pi f t) + Vdc

    Power. The power delivered through the electrode, V times the total
    current (conduction and displacement) from the electrode into the
    plasma, is averaged over updateCycles periods, and the amplitude is then
    changed at the start of a period (where the sine is zero):

        A <- A (power/measured power)^(relaxation/2),

    by at most the fraction maxChange. The exponent 1/2 corresponds to a
    power proportional to the square of the amplitude; the relaxation and
    the averaging allow for the plasma density following the amplitude over
    many periods.

    Blocking capacitor (blockingCapacitor):

        ideal     after every period the bias is changed by
                  -<I> T/capacitance, with <I> the mean current of that
                  period, until the mean current vanishes. Only the mean
                  current enters, so the capacitance is a gain here and can
                  be small for fast convergence;
        physical  Vdc = bias - Q/capacitance at every time step, with Q the
                  time integral of the current;
        none      Vdc = bias.

    \verbatim
    electrode
    {
        type            rfPowerElectrode;
        frequency       40e6;
        power           0.05;       // target mean power [W]
        amplitude       100;        // initial amplitude [V]
        updateCycles    5;          // periods per amplitude update
        relaxation      0.3;
        maxChange       0.1;
        blockingCapacitor ideal;
        capacitance     1e-12;      // [F]
        bias            0;          // initial (or fixed) bias [V]
        // area         1e-4;       // electrode area that the power refers
                                    // to [m2]; default: the patch area
        value           uniform 0;
    }
    \endverbatim

    The state (amplitude, bias, charge, running averages) is written with
    the field, so that a restarted run continues from it. One line per
    period is written to the log of the solver,

        rfPowerElectrode <patch>: t <time> A <amplitude> Vdc <bias>
            P <mean power> I <mean current>

    for the period that has ended.

    Needs the total current density field Jtot of somaFoam.

SourceFiles
    rfPowerElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef rfPowerElectrode_H
#define rfPowerElectrode_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rfPowerElectrode Declaration
\*---------------------------------------------------------------------------*/

class rfPowerElectrode
:
    public fixedValueFvPatchScalarField
{
    // Private data

        // Settings

            scalar frequency_;

            //- Target mean power [W]
            scalar power_;

            label updateCycles_;
            scalar relaxation_;
            scalar maxChange_;

            //- ideal, physical or none
            word blockingCapacitor_;
            scalar capacitance_;

            //- Electrode area that the power refers to; patch area if <= 0
            scalar area_;

        // State

            //- Amplitude [V] and bias [V]
            scalar amplitude_;
            scalar bias_;

            //- Capacitor charge [C] (physical capacitor)
            scalar charge_;

            //- Voltage applied in the last time step
            scalar lastVoltage_;

            //- Period in progress
            label cycle_;

            //- Sums over the period in progress: time, I dt
            scalar cycleTime_;
            scalar cycleCharge_;
            scalar cycleEnergy_;

            //- Sums since the last amplitude update: periods, time, V I dt
            label updatePeriods_;
            scalar updateTime_;
            scalar updateEnergy_;

            //- Time index of the last update of the sums
            label curTimeIndex_;


    // Private member functions

        //- Add the last time step to the sums; at the end of a period
        //  update the bias and, every updateCycles periods, the amplitude
        void advance();


public:

    //- Runtime type information
    TypeName("rfPowerElectrode");


    // Constructors

        //- Construct from patch and internal field
        rfPowerElectrode
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rfPowerElectrode
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rfPowerElectrode
        (
            const rfPowerElectrode&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rfPowerElectrode(const rfPowerElectrode&);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(new rfPowerElectrode(*this));
        }

        //- Construct as copy setting internal field reference
        rfPowerElectrode
        (
            const rfPowerElectrode&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>(new rfPowerElectrode(*this, iF));
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
