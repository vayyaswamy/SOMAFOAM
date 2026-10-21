/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
Foam::plasmaPotential

Description
    Prescribed electrode potential. The model sets the waveform:

        directCurrent           V = amplitude
        sinFrequencyModulated   V = amplitude sin(2 pi f t) + bias
        cosFrequencyModulated   V = amplitude cos(2 pi f t) + bias
        chirp                   V = amplitude e(t) sin(phi(t)) + bias
        table                   V = amplitude g(t) + bias

    chirp: the frequency f(t) is interpolated linearly in frequencyTable
    ((t f) ...) and held constant outside it; the phase is its exact integral,
    phi(t) = phase + 2 pi int_0^t f dt, so the waveform stays continuous when
    the frequency changes. The optional envelopeTable ((t e) ...) modulates
    the amplitude (default e = 1). phase is in degrees (default 0).

    table: g(t) is interpolated linearly in voltageTable ((t g) ...) or in
    the two-column text file voltageTableFile (t g per line), and held
    constant outside it. Use it for any waveform (triangular, pulsed,
    multi-harmonic, ...) generated elsewhere.

    Example:
        electrode
        {
            type            plasmaPotential;
            model           chirp;
            amplitude       uniform 100;
            frequencyTable  ((0 13.56e6) (50e-6 100e6));
            // envelopeTable ((0 1) (50e-6 0.5));
            bias            0;
            value           uniform 0;
        }

SourceFiles
plasmaPotential.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaPotential_H
#define plasmaPotential_H

#include "fixedValueFvPatchFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaPotential
:
public fixedValueFvPatchScalarField
{
    // Private data
    scalarField amplitude_;
    word modelName_;
    scalar frequency_;
    scalar bias_;

    //- chirp: frequency [Hz] against time [s], envelope, initial phase [deg]
    List<Tuple2<scalar, scalar> > frequencyTable_;
    List<Tuple2<scalar, scalar> > envelopeTable_;
    scalar phase_;

    //- table: waveform against time [s], given inline or read from a file
    List<Tuple2<scalar, scalar> > voltageTable_;
    fileName voltageTableFile_;


    // Private member functions

    //- Read the model-specific entries and check them
    void readModel(const dictionary& dict);

    //- Linear interpolation in a table sorted in time, held constant outside
    static scalar interpolate
    (
        const List<Tuple2<scalar, scalar> >& table,
        const scalar t
    );

    //- Phase of the chirp at time t [rad]
    scalar chirpPhase(const scalar t) const;

  public:

    //- Runtime type information
    TypeName("plasmaPotential");


    // Constructors

    //- Construct from patch and internal field
    plasmaPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    plasmaPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  plasmaPotential onto a new patch
    plasmaPotential
    (
     const plasmaPotential&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    plasmaPotential
    (
        const plasmaPotential&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
         new plasmaPotential(*this)
        );
    }

    //- Construct as copy setting internal field reference
    plasmaPotential
    (
         const plasmaPotential&,
         const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new plasmaPotential(*this, iF)
        );
    }


    // Member functions

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMap
    (
        const fvPatchFieldMapper&
    );

    //- Reverse map the given fvPatchField onto this fvPatchField
    virtual void rmap
    (
        const fvPatchScalarField&,
        const labelList&
    );

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
