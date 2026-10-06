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

        V(t) = A g(t) + Vdc

    Waveform g (waveform), periodic with the frequency f:

        sine       sin(2 pi f t)
        harmonics  sum over the entries (n a phi) of harmonics of
                   a sin(2 pi n f t + phi), phi in degrees
        table      linear interpolation in waveformTable ((x g) ...), with
                   x = f t modulo 1 the phase in the period, from 0 to 1

    Power. The target is constant (power) or follows powerTable
    ((t P) ...), interpolated linearly in time and evaluated at the middle
    of each period; with powerRepeat the table, which must start at t = 0,
    is repeated. The power delivered through the electrode, V times the
    total current (conduction and displacement) from the electrode into the
    plasma, is averaged over each period. No voltage is applied in periods
    with zero target. The amplitude A is set at the start of a period
    (tracking):

        feedback  every updateCycles periods
                      A <- A (target/measured)^(relaxation/2),
                  by at most the fraction maxChange; when the target
                  changes from one period to the next A is also scaled by
                  the square root of the ratio (feed-forward);
        learning  for a repeated powerTable: an amplitude is kept for each
                  period of the table and corrected as above every time
                  that period has been run, so that after some repetitions
                  the power follows the table through its fast changes.

    The exponent 1/2 corresponds to a power proportional to the square of
    the amplitude. A is kept below amplitudeMax.

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
        power           0.05;       // target mean power [W], or:
        // powerTable   ((0 0.05) (2.5e-6 0.05) (2.5e-6 0.01) (5e-6 0.01));
        // powerRepeat  yes;
        // tracking     learning;   // default feedback
        amplitude       100;        // initial amplitude [V]
        // amplitudeMax 1000;
        // waveform     harmonics;  // default sine
        // harmonics    ((1 1 0) (2 0.5 90));
        // waveform     table;
        // waveformTable ((0 0) (0.25 1) (0.75 -1) (1 0));
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

    The state (amplitude, bias, charge, running averages, learned
    amplitudes) is written with the field, so that a restarted run
    continues from it. One line per period is written to the log of the
    solver,

        rfPowerElectrode <patch>: t <time> A <amplitude> Vdc <bias>
            P <mean power> I <mean current> target <target power>

    for the period that has ended.

    Needs the total current density field Jtot of somaFoam.

SourceFiles
    rfPowerElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef rfPowerElectrode_H
#define rfPowerElectrode_H

#include "fixedValueFvPatchFields.H"
#include "Tuple2.H"
#include "vector.H"

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

            //- Target mean power [W]: constant, or a table in time
            scalar power_;
            List<Tuple2<scalar, scalar> > powerTable_;
            bool powerRepeat_;

            //- feedback or learning
            word tracking_;

            //- Waveform: sine, harmonics or table
            word waveform_;
            List<vector> harmonics_;
            List<Tuple2<scalar, scalar> > waveformTable_;

            scalar amplitudeMax_;

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
            //  and target power times time
            label updatePeriods_;
            scalar updateTime_;
            scalar updateEnergy_;
            scalar updateTarget_;

            //- Target power of the period in progress
            scalar currentTarget_;

            //- Amplitudes of the periods of a repeated power table
            //  (tracking learning)
            scalarField learned_;

            //- Time index of the last update of the sums
            label curTimeIndex_;


    // Private member functions

        //- Add the last time step to the sums; at the end of a period
        //  update the bias and, every updateCycles periods, the amplitude
        void advance();

        //- Waveform at a time
        scalar shape(const scalar t) const;

        //- Target power at a time
        scalar targetPower(const scalar t) const;

        //- Linear interpolation in a table of pairs
        static scalar interpolate
        (
            const List<Tuple2<scalar, scalar> >& table,
            const scalar x
        );

        //- Factor for the amplitude from a target and a measured power
        scalar correction(const scalar target, const scalar measured) const;


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
