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
    Foam::driftDiffusionNegativeIonDensity

Description
    Absorbing wall for a negative ion species in the drift-diffusion
    approximation, the mirror of driftDiffusionPositiveIonDensity. Where the
    ions drift to the wall (E.n < 0) the face value follows the cell value, so
    the drift carries them out; where they drift away from it the face value
    tends to zero and the wall does not inject them.

    A fixedValue 0 condition is not suitable with the Scharfetter-Gummel flux:
    the drift flux through a wall face is the drift velocity times the face
    value, so it blocks the drift out of the domain, and the ions pile up in
    the wall cells.

    Usage:
        ground
        {
            type    driftDiffusionNegativeIonDensity;
            value   uniform 0;
        }

SourceFiles
    driftDiffusionNegativeIonDensity.C

\*---------------------------------------------------------------------------*/

#ifndef driftDiffusionNegativeIonDensity_H
#define driftDiffusionNegativeIonDensity_H

#include "fvPatchFields.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class driftDiffusionNegativeIonDensity Declaration
\*---------------------------------------------------------------------------*/


class driftDiffusionNegativeIonDensity
:
    public mixedFvPatchField<scalar>
{
    // Private data
    scalar seec_;

    word fieldName_;

    public:

    //- Runtime type information
    TypeName("driftDiffusionNegativeIonDensity");


    // Constructors

        //- Construct from patch and internal field
        driftDiffusionNegativeIonDensity
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        driftDiffusionNegativeIonDensity
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inletOutletFvPatchField onto a new patch
        driftDiffusionNegativeIonDensity
        (
            const driftDiffusionNegativeIonDensity&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        driftDiffusionNegativeIonDensity
        (
            const driftDiffusionNegativeIonDensity&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new driftDiffusionNegativeIonDensity(*this)
            );
        }

        //- Construct as copy setting internal field reference
        driftDiffusionNegativeIonDensity
        (
            const driftDiffusionNegativeIonDensity&,
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
                new driftDiffusionNegativeIonDensity(*this, iF)
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
