/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::electronTemperatureWallFlux

Description
    Wall condition for the electron temperature in flux form, for use with
    the efullImplicit electron temperature model. The temperature itself has
    a zero gradient at the wall; the energy equation takes the total electron
    energy flux through the wall face from this condition:

        outgoing:  w k Te times the flux of plasma electrons to the wall,
                   w = 2 where they leave by their thermal motion and
                   w = 5/2 where they drift to the wall (Edepend and the
                   normal field pointing out of the wall), as in
                   driftDiffusionElectronDensity;
        incoming:  2 k Tse per secondary electron and 2 k TFN per
                   field-emitted electron.

    The loss is implicit in the temperature of the cell next to the wall
    with a coefficient that is never negative, so the condition is stable
    for any size of that cell. The electronTemperature condition expresses
    the same balance as a temperature gradient at the wall, which amplifies
    the wall temperature and fails when the cell is too coarse.

    The entries are those of electronTemperature:

    \verbatim
    electrode
    {
        type            electronTemperatureWallFlux;
        seec            0.05;
        Tse             1000;
        Edepend         true;
        TFN             1000;
        field_emission  false;
        field_enhancement_factor 55;
        work_function   5;
    }
    \endverbatim

SourceFiles
    electronTemperatureWallFlux.C

\*---------------------------------------------------------------------------*/

#ifndef electronTemperatureWallFlux_H
#define electronTemperatureWallFlux_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class electronTemperatureWallFlux Declaration
\*---------------------------------------------------------------------------*/

class electronTemperatureWallFlux
:
    public zeroGradientFvPatchScalarField
{
    // Private data

        //- Secondary emission coefficient and temperature of the secondary
        //  electrons [K]
        scalar seec_;
        scalar Tse_;

        //- Thermal flux to the wall only where the field repels the
        //  electrons from it
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
    TypeName("electronTemperatureWallFlux");


    // Constructors

        //- Construct from patch and internal field
        electronTemperatureWallFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronTemperatureWallFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        electronTemperatureWallFlux
        (
            const electronTemperatureWallFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronTemperatureWallFlux(const electronTemperatureWallFlux&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new electronTemperatureWallFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        electronTemperatureWallFlux
        (
            const electronTemperatureWallFlux&,
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
                new electronTemperatureWallFlux(*this, iF)
            );
        }


    // Member functions

        //- Flux of electrons emitted by the wall [1/m2/s]
        tmp<scalarField> emittedFlux() const;

        //- Energy flux brought in by the emitted electrons [W/m2]
        tmp<scalarField> emittedEnergyFlux() const;

        //- Energy carried to the wall per plasma electron, divided by its
        //  temperature [J/K]
        tmp<scalarField> energyPerElectron() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
