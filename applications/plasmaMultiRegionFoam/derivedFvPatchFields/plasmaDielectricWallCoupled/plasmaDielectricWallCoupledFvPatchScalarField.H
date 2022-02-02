/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    plasmaDielectricWallCoupledFvPatchScalarField

Description

    runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

SourceFiles
    plasmaDielectricWallCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaDielectricWallCoupledFvPatchScalarField_H
#define plasmaDielectricWallCoupledFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaDielectricWallCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

        const word neighbourPatchName_;
public:

    //- Runtime type information
    TypeName("plasmaDielectricWall");


    // Constructors

        //- Construct from patch and internal field
        plasmaDielectricWallCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plasmaDielectricWallCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        plasmaDielectricWallCoupledFvPatchScalarField
        (
            const plasmaDielectricWallCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaDielectricWallCoupledFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        plasmaDielectricWallCoupledFvPatchScalarField
        (
            const plasmaDielectricWallCoupledFvPatchScalarField&,
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
                new plasmaDielectricWallCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        const fvPatchScalarField& sigma() const;

/*        const fvPatchVectorField& E() const;*/

        const fvPatchScalarField& epsilon() const;

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
