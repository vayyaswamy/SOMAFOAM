/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::driftDiffusionElectronDensityFvPatchField

Description
    Foam::driftDiffusionElectronDensityFvPatchField

SourceFiles
    driftDiffusionElectronDensityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef driftDiffusionOpenBoundary_H
#define driftDiffusionOpenBoundary_H

#include "fvPatchFields.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class driftDiffusionElectronDensityFvPatch Declaration
\*---------------------------------------------------------------------------*/


class driftDiffusionOpenBoundary
:
    public mixedFvPatchField<scalar>
{
    // Private data
    scalar chargeSign_;

    word fieldName_;

    public:

    //- Runtime type information
    TypeName("driftDiffusionOpenBoundary");


    // Constructors

        //- Construct from patch and internal field
        driftDiffusionOpenBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        driftDiffusionOpenBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inletOutletFvPatchField onto a new patch
        driftDiffusionOpenBoundary
        (
            const driftDiffusionOpenBoundary&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        driftDiffusionOpenBoundary
        (
            const driftDiffusionOpenBoundary&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new driftDiffusionOpenBoundary(*this)
            );
        }

        //- Construct as copy setting internal field reference
        driftDiffusionOpenBoundary
        (
            const driftDiffusionOpenBoundary&,
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
                new driftDiffusionOpenBoundary(*this, iF)
            );
        }


    // Member functions

        

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<scalar>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "inletOutletFvPatchField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
