/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    plasmaDielectricEpsilonSlaveFvPatchScalarField

SourceFiles
    plasmaDielectricEpsilonSlaveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaDielectricEpsilonSlaveFvPatchScalarField_H
#define plasmaDielectricEpsilonSlaveFvPatchScalarField_H

#include "plasmaDielectricRegionCoupleBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaDielectricEpsilonFvPatchScalarField;
class coupledPotentialFvPatchScalarField;

/*---------------------------------------------------------------------------*\
        Class plasmaDielectricEpsilonSlaveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class plasmaDielectricEpsilonSlaveFvPatchScalarField
:
    public plasmaDielectricRegionCoupleBase
{
public:

    //- Runtime type information
    TypeName("plasmaDielectricEpsilonSlave");


    // Constructors

        //- Construct from patch and internal field
        plasmaDielectricEpsilonSlaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plasmaDielectricEpsilonSlaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  plasmaDielectricEpsilonSlaveFvPatchScalarField onto a new patch
        plasmaDielectricEpsilonSlaveFvPatchScalarField
        (
            const plasmaDielectricEpsilonSlaveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        plasmaDielectricEpsilonSlaveFvPatchScalarField
        (
            const plasmaDielectricEpsilonSlaveFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaDielectricEpsilonSlaveFvPatchScalarField(*this)
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
                new plasmaDielectricEpsilonSlaveFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return shadow patch field
            virtual const plasmaDielectricRegionCoupleBase& shadowPatchField() const;


        // Evaluation functions

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
