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

#include "momentumElectronDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumElectronDensity::momentumElectronDensity
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


Foam::momentumElectronDensity::momentumElectronDensity
(
    const momentumElectronDensity& ptf,
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


Foam::momentumElectronDensity::momentumElectronDensity
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



Foam::momentumElectronDensity::momentumElectronDensity
(
    const momentumElectronDensity& ptf
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



Foam::momentumElectronDensity::momentumElectronDensity
(
    const momentumElectronDensity& ptf,
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


void Foam::momentumElectronDensity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField n = patch().nf();

    const fvPatchField<scalar>& Uef=
        patch().lookupPatchField<volScalarField, scalar>("U_electron"); 

    c
    
 
 
    this->refValue() = 0.0;

    this->valueFraction() =  ;

    this->refGrad() = 0.0;

    mixedFvPatchField<scalar>::updateCoeffs();

   
}


void Foam::momentumElectronDensity::write(Ostream& os) const
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


void Foam::momentumElectronDensity::operator=
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
        momentumElectronDensity
    );
} // End namespace Foam

// ************************************************************************* //
