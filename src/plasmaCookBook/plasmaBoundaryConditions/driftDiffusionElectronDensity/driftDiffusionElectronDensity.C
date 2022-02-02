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

#include "driftDiffusionElectronDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftDiffusionElectronDensity::driftDiffusionElectronDensity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0),
    Edepend_(true),
    FE_(false),
    beta_(1.0),
    wf_(1.0)
    
{
    //Info << "Constructor 1" << endl;
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::driftDiffusionElectronDensity::driftDiffusionElectronDensity
(
    const driftDiffusionElectronDensity& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{
    //Info << "Constructor 3" << endl;
    //fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::driftDiffusionElectronDensity::driftDiffusionElectronDensity
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Edepend_(readBool(dict.lookup("Edepend"))),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
{
    //Info << "Constructor 2 " << endl;
    //Info << "seec_ " << seec_ << endl;
    //Info << "fieldName = " << iF.name() << endl;
    this->refValue() = 0.0;

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
    Info << "Edepend = " << Edepend_ << endl;
    Info << "FE_ = " << FE_ << endl;
    //Info << "work function = " << wf_ << endl;
    //Info << "field_enhancement_factor " << beta_ << endl; 
}



Foam::driftDiffusionElectronDensity::driftDiffusionElectronDensity
(
    const driftDiffusionElectronDensity& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{
    //Info << "Constructor 5" << endl;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::driftDiffusionElectronDensity::driftDiffusionElectronDensity
(
    const driftDiffusionElectronDensity& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{
    //Info << "Constructor 6" << endl;
    //Info << "seec = " << seec_ << endl;
    //fvPatchField<scalar>::operator=(this->patchInternalField());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::driftDiffusionElectronDensity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Info << "Inside updateCoeffs " << endl;
    vectorField n = patch().nf();

    const fvPatchField<scalar>& muef=
        patch().lookupPatchField<volScalarField, scalar>("mu_electron");

    const fvPatchField<scalar>& Def=
        patch().lookupPatchField<volScalarField, scalar>("D_electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te"); 

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    const fvPatchField<vector>& gradTef=
        patch().lookupPatchField<volVectorField, vector>("gradTe");

    const scalarField C1 = 0.25*sqrt(8.0*1.38e-23*Tef/9.1e-31/acos(-1.0)) + muef*(Ef&n) + 1.38e-23/1.602e-19*muef*(gradTef&n);

    const scalarField C2 = Def;

    scalarField Enorm = Ef&n ;

    scalarField Fifnorm = Fif&n;

    scalarField a = pos(mag(Enorm));   

    scalarField b = pos(Fifnorm);

    if (Edepend_)
    {
        a = pos(Enorm);    
    }

    const scalarField Gamma_se = seec_*(b*Fifnorm);

    scalarField c = pos(Enorm); // to ensure field emission happens only if E-field is pointing inward even if Edepend is set to false

    scalarField Gamma_FE = c*0.0;

    //Info << "FE_ = " << FE_ << endl;

    if (FE_)
    {
        //Info << "Enorm = " << Enorm << endl;

        scalarField vofy = 0.95 - sqr(3.79E-4)*beta_*c*mag(Enorm)*1E-2/sqr(wf_) ; // 1E-2 is for converting V/m to V/cm

        //Info << "vofy = " << vofy << endl;

        Gamma_FE = c*1.54E-6/1.602e-19*sqr(beta_*c*mag(Enorm) )/1.1/wf_*exp(-6.85E9*pow(wf_,1.5)*vofy/beta_/(c*mag(Enorm) + SMALL) ) ;

        //Info << "Gamma_FE = " << Gamma_FE << endl;
    }
    
 
 
    this->refValue() = 0.0;

    this->valueFraction() = a*C1/(C1+C2*this->patch().deltaCoeffs()) ;

    this->refGrad() = (Gamma_se + Gamma_FE)/(C2+SMALL);

    mixedFvPatchField<scalar>::updateCoeffs();

   
}


void Foam::driftDiffusionElectronDensity::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend")
        << Edepend_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_emission")
        << FE_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_enhancement_factor")
        << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("work_function")
        << wf_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::driftDiffusionElectronDensity::operator=
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
        driftDiffusionElectronDensity
    );
} // End namespace Foam

// ************************************************************************* //
