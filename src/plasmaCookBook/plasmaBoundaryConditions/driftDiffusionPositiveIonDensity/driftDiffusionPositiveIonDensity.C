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

#include "driftDiffusionPositiveIonDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftDiffusionPositiveIonDensity::driftDiffusionPositiveIonDensity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0.0),
    fieldName_(iF.name())
{
    //Info << "ion Constructor 1" << endl;
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::driftDiffusionPositiveIonDensity::driftDiffusionPositiveIonDensity
(
    const driftDiffusionPositiveIonDensity& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    fieldName_(ptf.fieldName_)
{}


Foam::driftDiffusionPositiveIonDensity::driftDiffusionPositiveIonDensity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    fieldName_(iF.name())
{
    //Info << "Constructor 2 " << endl;
    //Info << "fieldName = " << iF.name() << endl;
    this->refValue() = 0.0;

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::driftDiffusionPositiveIonDensity::driftDiffusionPositiveIonDensity
(
    const driftDiffusionPositiveIonDensity& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    fieldName_(ptf.fieldName_)
{}



Foam::driftDiffusionPositiveIonDensity::driftDiffusionPositiveIonDensity
(
    const driftDiffusionPositiveIonDensity& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    fieldName_(ptf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::driftDiffusionPositiveIonDensity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Info << "Inside updateCoeffs " << endl;
    vectorField n = patch().nf();

    

    //Info << "fieldName = " << fieldName_ << endl;

    const fvPatchField<scalar>& muf=
        patch().lookupPatchField<volScalarField, scalar>("mu_" + fieldName_);

    const fvPatchField<scalar>& Df=
        patch().lookupPatchField<volScalarField, scalar>("D_" + fieldName_);
       

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

   // const fvPatchField<vector>& Fif=
     //   patch().lookupPatchField<volVectorField, vector>("ionFlux");

    const scalarField C1 = -muf*(Ef&n);

    const scalarField C2 = Df;

    scalarField Enorm = Ef&n ;

    scalarField a = pos(Enorm);

    //Info << "Ef = " << Ef << endl;

    //Info << "n = " << n << endl;

    //Info << "Enorm = " << Enorm << endl;

    //Info << "a = " << a << endl;
 
    this->refValue() = 0.0;

    this->valueFraction() = (1-a)*C1/(C1+C2*this->patch().deltaCoeffs());

    //this->valueFraction() = 0.0;

    this->refGrad() = 0;

    mixedFvPatchField<scalar>::updateCoeffs();

    //Info << "Nef = " << Nef << endl;
}


void Foam::driftDiffusionPositiveIonDensity::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::driftDiffusionPositiveIonDensity::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );

    //Info << "operator " << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        driftDiffusionPositiveIonDensity
    );
} // End namespace Foam

// ************************************************************************* //
