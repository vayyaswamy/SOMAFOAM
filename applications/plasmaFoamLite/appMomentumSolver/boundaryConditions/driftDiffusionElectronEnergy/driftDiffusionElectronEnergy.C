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
#include "driftDiffusionElectronEnergy.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::driftDiffusionElectronEnergy::driftDiffusionElectronEnergy
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0.0),
    Tse_(0.0)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::driftDiffusionElectronEnergy::driftDiffusionElectronEnergy
(
    const driftDiffusionElectronEnergy& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
}

Foam::driftDiffusionElectronEnergy::driftDiffusionElectronEnergy
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Tse_(readScalar(dict.lookup("Tse")))
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::driftDiffusionElectronEnergy::driftDiffusionElectronEnergy
(
    const driftDiffusionElectronEnergy& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::driftDiffusionElectronEnergy::driftDiffusionElectronEnergy
(
    const driftDiffusionElectronEnergy& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::driftDiffusionElectronEnergy::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
  //Info << "Inside updateCoeffs " << endl;
    vectorField n = patch().nf();
    label patchi = this->patch().index();
    //Info << "fieldName = " << fieldName_ << endl;
    const fvPatchField<scalar>& muef=
        patch().lookupPatchField<volScalarField, scalar>("mobilityelectron");
    const fvPatchField<scalar>& nef=
        patch().lookupPatchField<volScalarField, scalar>("electron");
    const fvPatchField<scalar>& kef=
        patch().lookupPatchField<volScalarField, scalar>("kappa_electron");

    const fvPatchField<scalar>& Def=
        patch().lookupPatchField<volScalarField, scalar>("diffusionelectron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("eTemp");
    //const fvPatchField<vector>& gradTef=
    //    patch().lookupPatchField<volVectorField, scalar>("gradTe");    
    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("EField");
    const fvPatchField<vector>& vef=
        patch().lookupPatchField<volVectorField, vector>("velocityelectron");
    const volVectorField& vef_cell =
        db().objectRegistry::lookupObject<volVectorField>("velocityelectron");
    const volScalarField& nef_cell =
        db().objectRegistry::lookupObject<volScalarField>("electron");
    const volScalarField& tef_cell =
        db().objectRegistry::lookupObject<volScalarField>("eTemp");
    const volScalarField& kef_cell =
        db().objectRegistry::lookupObject<volScalarField>("kappa_electron");
    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    scalarField a = pos(vef_cell&n);

    const fvPatchVectorField& FiT = vef_cell.boundaryField()[patchi];
    scalarField temp1 = FiT&n;
    scalarField temp2 = FiT&n;
    scalarField temp3 = FiT&n;
    scalarField temp4 = FiT&n;
    forAll(a, facei)
    {
        //temp1[facei] = (vef_cell[facei]&n[facei])*(3/2)*nef_cell[facei]*1.38e-23 + 0.25*sqrt(8.0*1.38e-23*tef_cell[facei]/9.1e-31/acos(-1.0));//
        temp2[facei] = (vef[facei]&n[facei])*(3/2)*nef[facei]*1.38e-23 ;
        //temp3[facei] = - (3/5)*kef_cell[facei];
        //temp4[facei] = - (3/5)*kef[facei];
    }

    
    //const scalarField C2 = (5/3)*Def;
    //Info << "diffusion coefficient" << C2 << endl;
    //scalarField Enorm = Ef&n ;
    //scalarField Fifnorm = Fif&n;
    //scalarField a = pos();   
    //if (Edepend_)
    //{
        
    //}
    //const scalarField Gamma_se = seec_*(Fifnorm); 
 
    //this->refValue() = 0.0;
    //this->valueFraction() =0* (temp1+temp3*this->patch().deltaCoeffs())/(temp2+temp4*this->patch().deltaCoeffs()+1e-50);
    //this->refGrad() = 0;
    mixedFvPatchField<scalar>::updateCoeffs();
    
}


void Foam::driftDiffusionElectronEnergy::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::driftDiffusionElectronEnergy::operator=
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        driftDiffusionElectronEnergy
    );
} // End namespace Foam
// ************************************************************************* //
