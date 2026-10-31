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
    seec_(0.0),
    Edepend_(false)
   // fieldName_(iF.name())
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
    Edepend_(ptf.Edepend_)
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
    Edepend_(dict.lookup("Edepend"))
{
    //Info << "Constructor 2 " << endl;
    //Info << "fieldName = " << iF.name() << endl;
    this->refValue() = 0.0;

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::driftDiffusionElectronDensity::driftDiffusionElectronDensity
(
    const driftDiffusionElectronDensity& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_)
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
    Edepend_(ptf.Edepend_)
{
    //Info << "Constructor 6" << endl;
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

    

    //Info << "fieldName = " << fieldName_ << endl;

    const fvPatchField<scalar>& muef=
        patch().lookupPatchField<volScalarField, scalar>("mobilityelectron");

    const fvPatchField<scalar>& Def=
        patch().lookupPatchField<volScalarField, scalar>("diffusionelectron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("eTemp");

    //const fvPatchField<scalar>& Nef=
    //    patch().lookupPatchField<volScalarField, scalar>("N_electron");

    //const fvPatchField<vector>& gradTef=
    //    patch().lookupPatchField<volVectorField, scalar>("gradTe");    

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("EField");

    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    const scalarField C1 = 0.25*sqrt(8.0*1.38e-23*Tef/9.1e-31/acos(-1.0)) + muef*(Ef&n) ; //+ 1.38e-23/1.602e-19*muef*(gradTef&n);

    const scalarField C2 = Def;

    scalarField Enorm = Ef&n ;

    scalarField Fifnorm = Fif&n;

    scalarField a = pos(mag(Enorm));   

    if (Edepend_)
    {
        a = pos(Enorm);    
    }

    const scalarField Gamma_se = seec_*(Fifnorm);

    //Info << "Ef = " << Ef << endl;

    //Info << "n = " << n << endl;

    //Info << "Enorm = " << Enorm << endl;

    //Info << "a = " << a << endl;

    //Info << "b = " << b << endl; 
 
    this->refValue() = 0.0;

    this->valueFraction() = a*C1/(C1+C2*this->patch().deltaCoeffs()+1e-10*SMALL) ;

    this->refGrad() = a*Gamma_se/(C2+1e-10*SMALL);

    mixedFvPatchField<scalar>::updateCoeffs();

    //Info << "Nef = " << Nef << endl;
}


void Foam::driftDiffusionElectronDensity::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend")
        << Edepend_ << token::END_STATEMENT << nl;
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

    //Info << "operator " << endl;
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
