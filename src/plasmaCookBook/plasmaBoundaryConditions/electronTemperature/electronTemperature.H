/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::electronTemperature

Description
    Wall condition for the electron temperature, for use with the
    efullImplicit electron temperature model. The temperature itself has a
    zero gradient at the wall; the energy equation takes the total electron
    energy flux through the wall face from this condition:

        outgoing:  2 k Te per plasma electron that leaves by its thermal
                   motion (flux n vth/4) and 5/2 k Te per electron of the
                   flux beyond that, which drifts to the wall, as in
                   driftDiffusionElectronDensity;
        incoming:  2 k Tse per secondary electron and 2 k TFN per
                   field-emitted electron.

    The loss is implicit in the temperature of the cell next to the wall
    with a coefficient that is never negative, so the condition is stable
    for any size of that cell.

    \verbatim
    electrode
    {
        type            electronTemperature;
        seec            0.05;
        Tse             1000;
        TFN             1000;
        field_emission  false;
        field_enhancement_factor 55;
        work_function   5;
    }
    \endverbatim

    Edepend is accepted and not used.

SourceFiles
    electronTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef electronTemperature_H
#define electronTemperature_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class electronTemperature Declaration
\*---------------------------------------------------------------------------*/

class electronTemperature
:
    public zeroGradientFvPatchScalarField
{
    // Private data

        //- Secondary emission coefficient and temperature of the secondary
        //  electrons [K]
        scalar seec_;
        scalar Tse_;

        //- Not used (read for compatibility with electronTemperature)
        bool Edepend_;

        //- Temperature of the field-emitted electrons [K]
        scalar TFN_;

        //- Field emission, field enhancement factor and work function [eV]
        bool FE_;
        scalar beta_;
        scalar wf_;


    // Private member functions

        //- Secondary and field emission fluxes [1/m2/s]
        void emittedFluxes(scalarField& secondary, scalarField& field) const;


public:

    //- Runtime type information
    TypeName("electronTemperature");


    // Constructors

        //- Construct from patch and internal field
        electronTemperature
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronTemperature
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        electronTemperature
        (
            const electronTemperature&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronTemperature(const electronTemperature&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new electronTemperature(*this)
            );
        }

        //- Construct as copy setting internal field reference
        electronTemperature
        (
            const electronTemperature&,
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
                new electronTemperature(*this, iF)
            );
        }


    // Member functions

        //- Flux of electrons emitted by the wall [1/m2/s]
        tmp<scalarField> emittedFlux() const;

        //- Energy flux brought in by the emitted electrons [W/m2]
        tmp<scalarField> emittedEnergyFlux() const;

        //- Energy flux carried to the wall by a given flux of plasma
        //  electrons [1/m2/s], divided by their temperature [W/m2/K]:
        //  2 k for the thermal part of the flux, 5/2 k for the rest
        tmp<scalarField> energyFluxPerKelvin(const scalarField& flux) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
