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

#include "electronTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronTemperature::electronTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0),
    Tse_(0),
    Edepend_(false)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
}


Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_)
{}


Foam::electronTemperature::electronTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Tse_(readScalar(dict.lookup("Tse"))),
    Edepend_(dict.lookup("Edepend"))
{
    this->refValue() = 0;

    this->refGrad() = 0;
    this->valueFraction() = 0;

    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_)
{}



Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::electronTemperature::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField n = patch().nf();

    //Info << "electronTemperature updateCoeffs " << endl;

    const fvPatchField<scalar>& Nef=
        patch().lookupPatchField<volScalarField, scalar>("electron");

    const fvPatchField<scalar>& kappaef=
        patch().lookupPatchField<volScalarField, scalar>("eConductivity");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("eTemp");

    //Info << "Tef = " << Tef << endl;

    //Info << "Tecell = " << this->patchInternalField() << endl;

    //Info << "1 - f " << 1.0-this->valueFraction() << endl;

    //const fvPatchField<scalar>& gradTef=
    //    patch().lookupPatchField<volScalarField, scalar>("gradTe");    

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("EField");

    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

   // Info << "Step 1" << endl ;

    scalarField Enorm = Ef&n ;

    scalarField Fifnorm = Fif&n;

    scalarField a = pos(mag(Enorm));   

    if (Edepend_)
    {
        a = pos(Enorm);    
    }

    

    const scalarField C1 = 0.5*1.38e-23*0.25*Nef*sqrt(8.0*1.38e-23*Tef/9.1e-31/acos(-1.0)) - 2.5*1.38e-23*seec_*(Fif&n);

    const scalarField C2 = kappaef;

    const scalarField C3 = seec_*(Fifnorm)*1.38e-23*Tse_;

    this->refValue() = 0.0;

    //Info << "Step 2 " << endl;

    //Info << "C1 = " << C1 << endl;

    //Info << "C2 = "  << C2 << endl;

    //Info << "C3 = " << C3 << endl;

    this->valueFraction() = a*C1/(C1-C2*this->patch().deltaCoeffs()+1e-50) ;

    //Info << "value Fraction = " << this->valueFraction() << endl;

    this->refGrad() = a*C3/(C2+1e-50);

    //Info << "ref Grad = " << this->refGrad() << endl;

    

    mixedFvPatchField<scalar>::updateCoeffs();

    
}


void Foam::electronTemperature::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend")
        << Edepend_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::electronTemperature::operator=
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
        electronTemperature
    );
} // End namespace Foam

// ************************************************************************* //
