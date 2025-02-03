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

\*---------------------------------------------------------------------------*/

#include "dielectricSideWall.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dielectricSideWall::dielectricSideWall
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)    
{
    //Info << "Constructor 1" << endl;
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::dielectricSideWall::dielectricSideWall
(
    const dielectricSideWall& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    //Info << "Constructor 3" << endl;
    //fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::dielectricSideWall::dielectricSideWall
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    //Info << "Constructor 2 " << endl;
    //Info << "seec_ " << seec_ << endl;
    //Info << "fieldName = " << iF.name() << endl;
    this->refValue() = 0.0;

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::dielectricSideWall::dielectricSideWall
(
    const dielectricSideWall& ptf
)
:
    mixedFvPatchField<scalar>(ptf)
{
    //Info << "Constructor 5" << endl;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::dielectricSideWall::dielectricSideWall
(
    const dielectricSideWall& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    //Info << "Constructor 6" << endl;
    //Info << "seec = " << seec_ << endl;
    //fvPatchField<scalar>::operator=(this->patchInternalField());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dielectricSideWall::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Info << "Inside updateCoeffs " << endl;
    vectorField n = patch().nf();

    const fvPatchField<scalar>& surfC=
        patch().lookupPatchField<volScalarField, scalar>("surfC");
 
    this->refValue() = 0.0;

    this->valueFraction() = 0.0;

    this->refGrad() = 0.5*surfC/8.854e-12;

    mixedFvPatchField<scalar>::updateCoeffs();

   
}


void Foam::dielectricSideWall::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::dielectricSideWall::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dielectricSideWall
    );
} // End namespace Foam

// ************************************************************************* //
