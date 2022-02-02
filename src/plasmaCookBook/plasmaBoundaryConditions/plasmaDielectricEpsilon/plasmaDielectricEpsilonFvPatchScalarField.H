/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    plasmaDielectricEpsilonFvPatchScalarField

SourceFiles
    plasmaDielectricEpsilonFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaDielectricEpsilonFvPatchScalarField_H
#define plasmaDielectricEpsilonFvPatchScalarField_H

#include "plasmaDielectricRegionCoupleBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class plasmaDielectricEpsilonFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class plasmaDielectricEpsilonFvPatchScalarField
:
    public plasmaDielectricRegionCoupleBase
{
public:

    //- Runtime type information
    TypeName("plasmaDielectricEpsilon");


    // Constructors

        //- Construct from patch and internal field
        plasmaDielectricEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plasmaDielectricEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  plasmaDielectricEpsilonFvPatchScalarField onto a new patch
        plasmaDielectricEpsilonFvPatchScalarField
        (
            const plasmaDielectricEpsilonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        plasmaDielectricEpsilonFvPatchScalarField
        (
            const plasmaDielectricEpsilonFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaDielectricEpsilonFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaDielectricEpsilonFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            virtual void updateCoeffs();

            virtual void evaluate(const Pstream::commsTypes commsType);

            virtual tmp<scalarField> calcEpsilon
            (
                const plasmaDielectricRegionCoupleBase& owner,
                const plasmaDielectricRegionCoupleBase& neighbour
            ) const;

            virtual tmp<scalarField> calcPotential
            (
                const coupledPotentialFvPatchScalarField& owner,
                const coupledPotentialFvPatchScalarField& neighbour,
                const plasmaDielectricRegionCoupleBase& ownerEps
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
