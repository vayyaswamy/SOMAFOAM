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
#include "driftDiffusionElectronTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::driftDiffusionElectronTemperature::driftDiffusionElectronTemperature
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


Foam::driftDiffusionElectronTemperature::driftDiffusionElectronTemperature
(
    const driftDiffusionElectronTemperature& ptf,
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

Foam::driftDiffusionElectronTemperature::driftDiffusionElectronTemperature
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

Foam::driftDiffusionElectronTemperature::driftDiffusionElectronTemperature
(
    const driftDiffusionElectronTemperature& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::driftDiffusionElectronTemperature::driftDiffusionElectronTemperature
(
    const driftDiffusionElectronTemperature& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::driftDiffusionElectronTemperature::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
 vectorField n = patch().nf();//normal vector
    label patchi = this->patch().index();

    //getting values of different fields for the patch
    const volScalarField& muef =
        db().objectRegistry::lookupObject<volScalarField>("mobilityelectron");
    const volScalarField& Def =
        db().objectRegistry::lookupObject<volScalarField>("diffusionelectron");
    const volScalarField& electronCell =
        db().objectRegistry::lookupObject<volScalarField>("electron");
    const volScalarField& Tef =
        db().objectRegistry::lookupObject<volScalarField>("eTemp");
    //const volScalarField& Tef =
    //    db().objectRegistry::lookupObject<volScalarField>("eTemp");
    const volScalarField& cond =
        db().objectRegistry::lookupObject<volScalarField>("kappa_electron");
    const volVectorField& gradTe_cell =
        db().objectRegistry::lookupObject<volVectorField>("gradTe");
    const volVectorField& Fif =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux"); 
    const volVectorField& Ef =
        db().objectRegistry::lookupObject<volVectorField>("EField"); 

    const fvPatchField<vector>& vel_bound=
        patch().lookupPatchField<volVectorField, scalar>("velocityelectron");
    const fvPatchField<vector>& Ef_bound=
        patch().lookupPatchField<volVectorField, vector>("EField");
    const fvPatchField<vector>& fluxElec=
        patch().lookupPatchField<volVectorField, vector>("Flux_electron");
    const fvPatchField<scalar>& mob_bound=
        patch().lookupPatchField<volScalarField, scalar>("mobilityelectron");
    const fvPatchField<vector>& gradTe=
        patch().lookupPatchField<volVectorField, vector>("gradTe");
    const fvPatchField<scalar>& electron_n=
        patch().lookupPatchField<volScalarField, scalar>("electron");
         

    const fvPatchVectorField& FiT = Fif.boundaryField()[patchi];
    scalarField pp = FiT&n;
    scalarField qq = FiT&n;
    scalarField rr = FiT&n;
    forAll(pp, facei)
    {
    	label faceCelli = patch().faceCells()[facei];
        int a = pos(Ef[facei]&n[facei]);
        float vth = sqrt(8.0*1.38e-23*Tef[facei]/9.1e-31/acos(-1.0));
		pp[facei]= (5/2)*(a-1)*muef[facei]*(Ef[facei]&n[facei])*electronCell[facei]*1.38e-23 + (3/8)*vth*electronCell[facei]*1.38e-23;
        //qq[facei]=(3/2)*1.38e-23*Tse_*a*seec_*(Fif[facei]&n[facei])+ ((cond[facei]*gradTe_cell[facei] - cond[facei]*gradTe[facei]) & n[facei]);
        rr[facei]= (5/2)*1.38e-23*(fluxElec[facei]& n[facei]);
        Info <<"electron flux "<< fluxElec[facei]<<"pp" << pp[facei] <<"qq" << qq[facei] <<"rr" << rr[facei] << endl;
	}
    this->refValue() = 0;
    this->valueFraction() = 1-pp/(rr+1e-50); 
    this->refGrad() = 0*qq/(pp+1e-50)*this->patch().deltaCoeffs();
    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::driftDiffusionElectronTemperature::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::driftDiffusionElectronTemperature::operator=
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
        driftDiffusionElectronTemperature
    );
} // End namespace Foam
// ************************************************************************* //
