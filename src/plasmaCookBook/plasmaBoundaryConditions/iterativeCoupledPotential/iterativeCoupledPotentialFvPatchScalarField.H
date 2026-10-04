/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::iterativeCoupledPotentialFvPatchScalarField

Description
    Potential on a plasma-dielectric interface (regionCouple patch pair) for
    the segregated solver plasmaMultiRegionFoam, which solves the plasma and
    the dielectric regions in turn and iterates until the interface
    conditions are met (somaFoam solves the regions in one matrix and uses
    coupledPotential instead).

    One side of the interface fixes the interface potential and the other
    side imposes the normal gradient that satisfies Gauss's law,

        eps_d dPhi/dn_d + eps_p dPhi/dn_p = sigma

    with the outward normal of each region and the surface charge sigma. The
    solver compares the potential the gradient side obtains on the interface
    with the fixed one and relaxes the fixed one towards it.

    The permittivity is that of the cells next to the interface, taken from
    the field epsilonEff of the region if the solver provides it (the
    semi-implicit Poisson equation) and from epsilon otherwise.

    \verbatim
    dielectric_1                       // in 0/Phi (plasma)
    {
        type            iterativeCoupledPotential;
        remoteField     Phi;
        surfaceCharge   on;            // this region holds surfC
        fixesValue      yes;           // default: the master patch
        value           uniform 0;
    }
    plasma_1                           // in 0/dielectric_1/Phi
    {
        type            iterativeCoupledPotential;
        remoteField     Phi;
        surfaceCharge   off;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    iterativeCoupledPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef iterativeCoupledPotentialFvPatchScalarField_H
#define iterativeCoupledPotentialFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "regionCoupleFvPatch.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class iterativeCoupledPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class iterativeCoupledPotentialFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the potential in the other region
        word remoteFieldName_;

        //- The surface charge field surfC is on this side
        Switch surfaceCharge_;

        //- This side fixes the interface potential; the other side
        //  imposes the gradient
        Switch fixesValue_;


    // Private Member Functions

        //- Permittivity of the cells next to a patch
        static tmp<scalarField> permittivity(const fvPatch&);


public:

    //- Runtime type information
    TypeName("iterativeCoupledPotential");


    // Constructors

        //- Construct from patch and internal field
        iterativeCoupledPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        iterativeCoupledPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        iterativeCoupledPotentialFvPatchScalarField
        (
            const iterativeCoupledPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        iterativeCoupledPotentialFvPatchScalarField
        (
            const iterativeCoupledPotentialFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new iterativeCoupledPotentialFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        iterativeCoupledPotentialFvPatchScalarField
        (
            const iterativeCoupledPotentialFvPatchScalarField&,
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
                new iterativeCoupledPotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Does this side fix the interface potential?
        bool fixesValue() const
        {
            return fixesValue_;
        }

        //- Does this side hold the surface charge?
        bool surfaceCharge() const
        {
            return surfaceCharge_;
        }

        //- Surface charge of the interface, from whichever side holds it
        tmp<scalarField> sigma() const;

        //- The same patch field in the other region
        const iterativeCoupledPotentialFvPatchScalarField& shadow() const;

        //- Interface potential obtained by the other region, on this patch
        tmp<scalarField> shadowInterfaceValue() const;

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
