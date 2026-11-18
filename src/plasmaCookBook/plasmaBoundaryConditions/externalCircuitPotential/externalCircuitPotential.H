/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::externalCircuitPotential

Description
    Potential of an electrode connected to a voltage source through a
    resistor R and a capacitor C in series:

        Vs(t) = V + R I + Q_C/C,      dQ_C/dt = I,

    with V the potential of the electrode and I the current from the source
    into the electrode. The source is

        directCurrent           Vs = amplitude (1 - exp(-t/rampTime)) + bias
                                (no ramp with rampTime 0)
        sinFrequencyModulated   Vs = amplitude sin(2 pi frequency t) + bias
        cosFrequencyModulated   V  = amplitude cos(2 pi frequency t) + bias,
                                applied directly, without the circuit

    The current is the change of the charge Q on the electrode plus the
    conduction current Ic of the charged species from the electrode into the
    discharge,

        I = dQ/dt + Ic,       Q = sum over the faces of eps dPhi/dn A,

    so it includes the displacement current. The circuit is solved
    implicitly (backward Euler): with Q linearised about the present
    iterate, Q = Q* + Cg (V - V*), and Q^n, Q_C^n the values of the last
    time step,

        V (1 + Reff Cg/dt) = Vs - Q_C^n/C
                           - Reff [(Q* - Cg V* - Q^n)/dt + Ic],

    Reff = R + dt/C. The circuit is solved again only when the potential
    next to the electrode has changed since the last evaluation, so once
    per Poisson solve. The iterations of the Poisson equation within a step
    converge to the implicit solution, which is stable for any R C/dt; with
    one Poisson solve per step on a vacuum gap it is exactly the backward
    Euler step of the RC circuit. Cg only sets the convergence of the
    iterations and needs to be accurate within a factor of two; it is
    gapCapacitance if given, and otherwise, on a one-dimensional mesh, the
    capacitance of the plasma mesh between the electrode and the opposite
    end, A/sum(dx/eps).

    Before 2026-11 the current was taken from the previous step and the
    capacitor charge was never updated, so the displacement current of the
    first step was lost, a series capacitor had no effect, and a resolved
    displacement current would have made the update unstable for
    R Cg > dt/2.

    Example:
        electrode
        {
            type            externalCircuitPotential;
            model           directCurrent;
            amplitude       uniform 400;
            R               2e5;          // [Ohm]
            C               1e20;         // [F], optional, default: none
            rampTime        0;
            bias            0;
            // gapCapacitance 1.77e-13;   // [F], optional
            value           uniform 0;
        }

    The state (capacitor charge QC, electrode charge Qe, current I) is
    written with the field and read at a restart.

SourceFiles
    externalCircuitPotential.C

\*---------------------------------------------------------------------------*/

#ifndef externalCircuitPotential_H
#define externalCircuitPotential_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class externalCircuitPotential
:
    public fixedValueFvPatchScalarField
{
    // Private data

        scalarField amplitude_;
        word modelName_;
        scalar frequency_;
        scalar bias_;
        scalar R_;
        scalar C_;

        //- Time constant of the ramp of the source voltage (directCurrent)
        scalar rampTime_;

        //- Capacitance of the gap for the linearisation (<= 0: computed)
        scalar gapCapacitance_;

        //- Charge of the series capacitor at the start of the step [C]
        scalar QC_;

        //- Charge on the electrode at the start of the step [C]
        scalar Qe_;

        //- Current from the source at the last evaluation [A]
        scalar I_;

        //- Time index of the start-of-step values
        label timeIndex_;

        //- Has Qe_ been set (from the dictionary or a step)?
        bool QeSet_;

        //- Sum of the potentials of the cells next to the patch at the
        //  last evaluation; the circuit is solved again only after the
        //  Poisson equation has changed them, since the charge of the
        //  electrode is consistent with the potential of the patch only
        //  for the field solved with that potential
        scalar lastPhiP_;


    // Private member functions

        //- Charge on the electrode for the present field
        scalar electrodeCharge() const;

        //- Conduction current from the electrode into the discharge
        scalar conductionCurrent() const;

        //- Capacitance used to linearise the charge
        scalar gapCapacitance() const;

        //- Source voltage at time t
        scalar sourceVoltage(const scalar t) const;


public:

    //- Runtime type information
    TypeName("externalCircuitPotential");


    // Constructors

        externalCircuitPotential
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        externalCircuitPotential
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        externalCircuitPotential
        (
            const externalCircuitPotential&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        externalCircuitPotential(const externalCircuitPotential&);

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new externalCircuitPotential(*this)
            );
        }

        externalCircuitPotential
        (
            const externalCircuitPotential&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new externalCircuitPotential(*this, iF)
            );
        }


    // Member functions

        //- Current from the source into the electrode [A]
        scalar current() const
        {
            return I_;
        }

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
