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
    Edepend_(false),
    TFN_(0),
    FE_(false),
    beta_(1.0),
    wf_(1.0)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;

    //Info << "Constructor 1 for Te " << endl;
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
    Edepend_(ptf.Edepend_),
    TFN_(ptf.Tse_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
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
    Edepend_(readBool(dict.lookup("Edepend"))),
    TFN_(readScalar(dict.lookup("TFN"))),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
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
    Edepend_(ptf.Edepend_),
    TFN_(ptf.TFN_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
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
    Edepend_(ptf.Edepend_),
    TFN_(ptf.TFN_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
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
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    const fvPatchField<scalar>& kappaef=
        patch().lookupPatchField<volScalarField, scalar>("kappa_electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te");

    //Info << "Tef = " << Tef << endl;

    //Info << "kappef = " << kappaef << endl;

    //Info << "Nef = " << Nef << endl;

    //Info << "Tecell = " << this->patchInternalField() << endl;

    //Info << "1 - f " << 1.0-this->valueFraction() << endl;

    //const fvPatchField<scalar>& gradTef=
    //    patch().lookupPatchField<volScalarField, scalar>("gradTe");    

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

   // Info << "Step 1" << endl ;

    scalarField Enorm = Ef&n ;

    scalarField Fifnorm = Fif&n;

    scalarField a = pos(mag(Enorm));  

    scalarField b = pos(Fifnorm);

    if (Edepend_)
    {
        a = pos(Enorm);    
    }

    scalarField c = pos(Enorm); // to ensure field emission happens only if E-field is pointing inward even if Edepend is set to false

    const scalarField Gamma_se = seec_*(b*Fifnorm);

    scalarField Gamma_FE = 0.0*c;

    if (FE_)
    {
        //Info << "Enorm = " << Enorm << endl;

        scalarField vofy = 0.95 - sqr(3.79E-4)*beta_*c*mag(Enorm)*1E-2/sqr(wf_) ; // 1E-2 is for converting V/m to V/cm

        //Info << "vofy = " << vofy << endl;

        Gamma_FE = 1.54E-6/1.602e-19*sqr(beta_*c*mag(Enorm) )/1.1/wf_*exp(-6.85E9*pow(wf_,1.5)*vofy/beta_/(c*mag(Enorm) + SMALL) ) ;

        //Info << "Gamma_FE = " << Gamma_FE << endl;
    }
    

    const scalarField C1 = 0.5*1.38e-23*0.25*Nef*sqrt(8.0*1.38e-23*Tef/9.1e-31/acos(-1.0)) - 2.5*1.38e-23*(Gamma_se + Gamma_FE);

    const scalarField C2 = kappaef;

    const scalarField C3 = 1.38e-23*(Gamma_se*Tse_ + Gamma_FE*TFN_);

    this->refValue() = 0.0;

    

    this->valueFraction() = a*C1/(C1-C2*this->patch().deltaCoeffs()) ;

    //Info << "value Fraction = " << this->valueFraction() << endl;

    this->refGrad() = C3/(C2+SMALL);

    //Info << "ref Grad = " << this->refGrad() << endl;

    //Info << "value fraction = " << this->valueFraction() << endl;
    

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
     os.writeKeyword("TFN")
        << TFN_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_emission")
        << FE_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_enhancement_factor")
        << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("work_function")
        << wf_ << token::END_STATEMENT << nl;
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
