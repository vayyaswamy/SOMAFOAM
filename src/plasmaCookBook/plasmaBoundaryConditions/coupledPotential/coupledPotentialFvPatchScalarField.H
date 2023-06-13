/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    coupledPotentialFvPatchScalarField

SourceFiles
    coupledPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPotentialFvPatchScalarField_H
#define coupledPotentialFvPatchScalarField_H

#include "regionCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class coupledPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledPotentialFvPatchScalarField
:
    public regionCouplingFvPatchScalarField
{
    // Private data

        Switch surfaceCharge_;

public:

    //- Runtime type information
    TypeName("coupledPotential");


    // Constructors

        //- Construct from patch and internal field
        coupledPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given coupledPotentialFvPatchScalarField
        //  onto a new patch
        coupledPotentialFvPatchScalarField
        (
            const coupledPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coupledPotentialFvPatchScalarField
        (
            const coupledPotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coupledPotentialFvPatchScalarField(*this)
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
                new coupledPotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            const Switch& surfaceCharge() const
            {
                return surfaceCharge_;
            }

            Switch& surfaceCharge()
            {
                return surfaceCharge_;
            }

            //- Return shadow patch field
            virtual const coupledPotentialFvPatchScalarField&
                shadowPatchField() const;

            //- Return neighbour field given internal cell data
            virtual tmp<scalarField> patchNeighbourField() const;

            //- Return wall potential
            virtual tmp<scalarField> Phiw() const;

            //- Return face-cell potential
            virtual tmp<scalarField> Phic() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Calculate source term
            virtual tmp<scalarField> source() const;

            //- Manipulate matrix
            virtual void manipulateMatrix(fvScalarMatrix& matrix);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
