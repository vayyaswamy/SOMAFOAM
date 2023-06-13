/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    plasmaDielectricRegionCoupleBase

SourceFiles
    plasmaDielectricRegionCoupleBase.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaDielectricRegionCoupleBase_H
#define plasmaDielectricRegionCoupleBase_H

#include "regionCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class coupledPotentialFvPatchScalarField;

/*---------------------------------------------------------------------------*\
                  Class plasmaDielectricRegionCoupleBase Declaration
\*---------------------------------------------------------------------------*/

class plasmaDielectricRegionCoupleBase
:
    public regionCouplingFvPatchScalarField
{
    // Private data


public:

    //- Runtime type information
    TypeName("plasmaDielectricRegionCoupleBase");


    // Constructors

        //- Construct from patch and internal field
        plasmaDielectricRegionCoupleBase
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plasmaDielectricRegionCoupleBase
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given plasmaDielectricRegionCoupleBase onto a new patch
        plasmaDielectricRegionCoupleBase
        (
            const plasmaDielectricRegionCoupleBase&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        plasmaDielectricRegionCoupleBase
        (
            const plasmaDielectricRegionCoupleBase&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaDielectricRegionCoupleBase(*this)
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
                new plasmaDielectricRegionCoupleBase(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return false to disable interpolation on this patch
            virtual bool coupled() const
            {
                return false;
            }

            //- Return shadow patch field
            const plasmaDielectricRegionCoupleBase& shadowPatchField() const;

            //- Return original patch field
            const scalarField& originalPatchField() const
            {
                return regionCouplingFvPatchScalarField::originalPatchField();
            }

            const scalarField& originalPatchInternalField() const
            {
                return this->patchInternalField();
            }

            virtual tmp<scalarField> forig() const;

            virtual tmp<scalarField> korig() const;

            virtual tmp<scalarField> kc() const;

            virtual tmp<scalarField> kw() const;



        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

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

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
