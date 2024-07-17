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
#include "electronDensity_wilsonBC3.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronDensity_wilsonBC3::electronDensity_wilsonBC3
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0.0)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::electronDensity_wilsonBC3::electronDensity_wilsonBC3
(
    const electronDensity_wilsonBC3& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_)
{
}


Foam::electronDensity_wilsonBC3::electronDensity_wilsonBC3
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec")))
{
    this->refValue() = 0.0;

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::electronDensity_wilsonBC3::electronDensity_wilsonBC3
(
    const electronDensity_wilsonBC3& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_)
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::electronDensity_wilsonBC3::electronDensity_wilsonBC3
(
    const electronDensity_wilsonBC3& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::electronDensity_wilsonBC3::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    vectorField n = patch().nf();//normal vector
    label patchi = this->patch().index();

    //getting values of different fields for the patch
    //const volScalarField& muef =
    //    db().objectRegistry::lookupObject<volScalarField>("mobilityelectron");
    //const volScalarField& Def =
    //    db().objectRegistry::lookupObject<volScalarField>("diffusionelectron");
    const volScalarField& Tef =
        db().objectRegistry::lookupObject<volScalarField>("eTemp");
    const volScalarField& electrn =
        db().objectRegistry::lookupObject<volScalarField>("electron");
    const volVectorField& Ef =
        db().objectRegistry::lookupObject<volVectorField>("EField");
    const volVectorField& Fif =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux"); 

    const fvPatchField<vector>& vel_bound=
        patch().lookupPatchField<volVectorField, scalar>("velocityelectron");
    const fvPatchField<vector>& Ef_bound=
        patch().lookupPatchField<volVectorField, scalar>("EField");
    const fvPatchField<scalar>& mob_bound=
        patch().lookupPatchField<volScalarField, scalar>("mobilityelectron");
    const fvPatchField<scalar>& muef=
        patch().lookupPatchField<volScalarField, scalar>("mobilityelectron");
    const fvPatchField<scalar>& Def=
        patch().lookupPatchField<volScalarField, scalar>("diffusionelectron");


    const fvPatchVectorField& FiT = Fif.boundaryField()[patchi];
    //scalarField expterm = FiT&n;
    scalarField pp = FiT&n;
    scalarField qq = FiT&n;
    forAll(pp, facei)
    {
    	label faceCelli = patch().faceCells()[facei];
        int a = pos(Ef[faceCelli]&n[facei]);
        float vth = sqrt(8.0*1.38e-23*Tef[faceCelli]/9.1e-31/acos(-1.0)) ;
		//expterm[facei]= a*seec_*(Fif[faceCelli]&n[facei]);//0.5*vth*(1-a)*seec_*(Fif[facei]&n[facei])/(muef[facei]*(Ef_bound[facei]&n[facei]))+2*(1-a)*seec_*(Fif[facei]&n[facei]);
        pp[facei] = (a-1)*muef[facei]*(Ef[faceCelli]&n[facei]) + 0.5*vth - Def[facei]*this->patch().deltaCoeffs()[facei];//(1-2*a)*muef[facei]*(Ef[facei]&n[facei]) + 0.5*vth- Def[facei]*this->patch().deltaCoeffs()[facei];
        qq[facei] = -muef[facei]*(Ef_bound[facei]&n[facei]) - Def[facei]*this->patch().deltaCoeffs()[facei];
        //int a = pos(-Ef[faceCelli]&n[facei]);
        //float vth = sqrt(8.0*1.38e-23*Tef[faceCelli]/9.1e-31/acos(-1.0)) ;
		//expterm[facei]= 0.5*vth*(1-a)*seec_*(Fif[faceCelli]&n[facei])/(muef[facei]*(Ef[faceCelli]&n[facei]))+2*(1-a)*seec_*(Fif[facei]&n[facei]);
        //pp[facei] = (1-2*a)*muef[facei]*(Ef[faceCelli]&n[facei]) + 0.5*vth- Def[facei]*this->patch().deltaCoeffs()[facei];
        //qq[facei] = -muef[facei]*(Ef_bound[facei]&n[facei]) - Def[facei]*this->patch().deltaCoeffs()[facei];
	}
    //Info << "lets see the value of muef" << muef->patchInternalField() << endl;
    this->valueFraction() = 1-pp/qq; 
    this->refValue() = 0;
    this->refGrad() = 0;//expterm/pp*this->patch().deltaCoeffs();
    mixedFvPatchField<scalar>::updateCoeffs();
}

void Foam::electronDensity_wilsonBC3::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void Foam::electronDensity_wilsonBC3::operator=
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
          - this->refGrad()
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electronDensity_wilsonBC3
    );
} // End namespace Foam
// ************************************************************************* //
