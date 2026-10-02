/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::electronMomentumTemperature

Description
    Electron temperature wall condition for electrons transported with the
    momentum model.

    Unlike electronTemperature, which assumes the drift-diffusion wall flux,
    the electron wall flux is taken from the solved electron velocity,
    Gamma = N_electron*(U_electron & n), so it stays consistent with
    whatever the velocity condition (e.g. electronThermalVelocity) imposes.

    The total wall energy flux is set to
        Q = eps k Te Gamma_out - k (Tse Gamma_se + TFN Gamma_FE)
    with Gamma_out = Gamma + Gamma_se + Gamma_FE the outgoing electron flux
    and eps the energy (in k Te) carried per outgoing electron:
        eps = 2                  where the field repels electrons from the
                                 wall (thermal loss)
        eps = driftEnergyFactor  where it drives them into the wall (drift
                                 loss); the default 2.5 matches
                                 electronTemperature, 2 gives a single
                                 Hagelaar-style closure
    With Edepend false the thermal branch is used everywhere. Edepend must
    match the electron velocity condition.
    The energy equation already advects 2.5 k Te Gamma through the wall, so
    the condition supplies the remaining conductive flux:
        kappa dTe/dn = k ((2.5 - eps) Gamma_out - 2.5 Gamma_inj) Te
                     + k (Tse Gamma_se + TFN Gamma_FE)

    Example:
    \verbatim
        electrode
        {
            type                     electronMomentumTemperature;
            seec                     0.05;
            Edepend                  true;   // optional, default true
            driftEnergyFactor        2.5;    // optional, default 2.5
            Tse                      1000;
            TFN                      1000;
            field_emission           false;
            field_enhancement_factor 55.0;
            work_function            5.0;
        }
    \endverbatim

SourceFiles
    electronMomentumTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef electronMomentumTemperature_H
#define electronMomentumTemperature_H

#include "fvPatchFields.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class electronMomentumTemperature Declaration
\*---------------------------------------------------------------------------*/

class electronMomentumTemperature
:
    public mixedFvPatchField<scalar>
{
    // Private data

        //- Secondary electron emission coefficient
        scalar seec_;

        //- Temperature of secondary electrons
        scalar Tse_;

        //- Temperature of field-emitted electrons
        scalar TFN_;

        //- Thermal loss only where the field repels electrons
        bool Edepend_;

        //- Energy (in k Te) per electron lost by drift
        scalar driftEnergyFactor_;

        //- Field emission switch and Fowler-Nordheim parameters
        bool FE_;
        scalar beta_;
        scalar wf_;


public:

    //- Runtime type information
    TypeName("electronMomentumTemperature");


    // Constructors

        //- Construct from patch and internal field
        electronMomentumTemperature
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronMomentumTemperature
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given electronMomentumTemperature
        //  onto a new patch
        electronMomentumTemperature
        (
            const electronMomentumTemperature&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronMomentumTemperature
        (
            const electronMomentumTemperature&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new electronMomentumTemperature(*this)
            );
        }

        //- Construct as copy setting internal field reference
        electronMomentumTemperature
        (
            const electronMomentumTemperature&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new electronMomentumTemperature(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<scalar>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
