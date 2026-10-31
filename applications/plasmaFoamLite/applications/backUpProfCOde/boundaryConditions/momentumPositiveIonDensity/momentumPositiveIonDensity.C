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

#include "boundaryConditions/momentumPositiveIonDensity/momentumPositiveIonDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::momentumPositiveIonDensity::momentumPositiveIonDensity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    fieldName_(iF.name())
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::momentumPositiveIonDensity::momentumPositiveIonDensity
(
    const momentumPositiveIonDensity& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_)
{}


Foam::momentumPositiveIonDensity::momentumPositiveIonDensity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    fieldName_(iF.name())
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::momentumPositiveIonDensity::momentumPositiveIonDensity
(
    const momentumPositiveIonDensity& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    fieldName_(ptf.fieldName_)
{}



Foam::momentumPositiveIonDensity::momentumPositiveIonDensity
(
    const momentumPositiveIonDensity& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    fieldName_(ptf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//============================================================================
void Foam::momentumPositiveIonDensity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField n = patch().nf();//normal vector
    label patchi = this->patch().index();

    //getting values of different fields for the patch
    const volVectorField& vef =
        db().objectRegistry::lookupObject<volVectorField>("velocity"+ fieldName_);
    const volVectorField& Ef =
        db().objectRegistry::lookupObject<volVectorField>("EField");
   
   
   
   
    const fvPatchVectorField& FiT = Ef.boundaryField()[patchi];
    scalarField temp1 = FiT&n;
    forAll(temp1, facei)
    {
        label faceCelli = patch().faceCells()[facei];
        int a = pos(mag(vef[facei]&n[facei]));
                temp1[facei]=a;//*vef[facei]&n[facei];
        }

    this->refValue() = 0;
    this->valueFraction() =1 - temp1;//;temp1/((vef_bound&n)+1e-8); 
    this->refGrad() = 0; 
    mixedFvPatchField<scalar>::updateCoeffs();
}

void Foam::momentumPositiveIonDensity::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void Foam::momentumPositiveIonDensity::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::operator=
    (
         this->valueFraction()*this->refValue()
         + (1.0 - this->valueFraction())*
        (
            this->patchInternalField()
            - this->refGrad()/this->patch().deltaCoeffs()
        )
    );
}
//-----------------------------------------------------------------------------//


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        momentumPositiveIonDensity
    );
} // End namespace Foam
// ************************************************************************* //
