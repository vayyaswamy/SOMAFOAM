/*---------------------------------------------------------------------------*\
Copyright (C) 2023 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
Foam::externalCircuit1D 

Description
    Class for electric potential boundary conditions.

SourceFiles
externalCircuit1D.C

\*---------------------------------------------------------------------------*/

#ifndef externalCircuit1D_H
#define externalCircuit1D_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class externalCircuit1D
    :
    public fixedValueFvPatchScalarField
    {
        private:

            scalarField amplitude_;
            word modelName_;
            scalar frequency_;
            scalar bias_;
            scalar R_;
            scalar distance_;
            bool init_;
            scalar oldValue_;
            scalar newValue_; 

        public:

            TypeName("externalCircuit1D");

            // Construct from patch and internal field
            externalCircuit1D
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&
            );

            // Construct from patch, internal field and dictionary
            externalCircuit1D
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const dictionary&            
            );

            // Construct by mapping given externalCircuit1D onto a 
            // new patch 
            externalCircuit1D
            (
                const externalCircuit1D&,
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const fvPatchFieldMapper&
            );

            // Construct as copy
            externalCircuit1D
            (
                const externalCircuit1D&
            );

            // Construct and return a clone
            virtual tmp<fvPatchScalarField> clone() const
            {
                return tmp<fvPatchScalarField>
                (
                    new externalCircuit1D(*this)
                );
            }

            // Construct a copy setting internal field reference
            externalCircuit1D
            (
                const externalCircuit1D&,
                const DimensionedField<scalar, volMesh>&
            );

            // Construct and return a clone setting internal field reference 
            virtual tmp<fvPatchScalarField> 
            clone(const DimensionedField<scalar, volMesh>& iF) const
            {
                return tmp<fvPatchScalarField>
                (
                    new externalCircuit1D(*this, iF)
                );
            }

            // Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            // Write @ I/O operations
            virtual void write(Ostream&) const;
    };
}

#endif