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
#include "electronTemperatureWilsonBC5.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::electronTemperatureWilsonBC5::electronTemperatureWilsonBC5
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


Foam::electronTemperatureWilsonBC5::electronTemperatureWilsonBC5
(
    const electronTemperatureWilsonBC5& ptf,
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

Foam::electronTemperatureWilsonBC5::electronTemperatureWilsonBC5
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

Foam::electronTemperatureWilsonBC5::electronTemperatureWilsonBC5
(
    const electronTemperatureWilsonBC5& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::electronTemperatureWilsonBC5::electronTemperatureWilsonBC5
(
    const electronTemperatureWilsonBC5& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::electronTemperatureWilsonBC5::updateCoeffs()
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
    //const volScalarField& kappae =
    //    db().objectRegistry::lookupObject<volScalarField>("diffusionelectron");
    const volScalarField& Tef =
        db().objectRegistry::lookupObject<volScalarField>("eTemp");
    const volScalarField& eDensit =
        db().objectRegistry::lookupObject<volScalarField>("eDensity");
    const volVectorField& Ef =
        db().objectRegistry::lookupObject<volVectorField>("EField");
    const volVectorField& Fif =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux"); 
    const volVectorField& flux_cell =
        db().objectRegistry::lookupObject<volVectorField>("Flux_electron"); 
    const fvPatchField<vector>& flux_bound=
        patch().lookupPatchField<volVectorField, scalar>("Flux_electron");
    const fvPatchField<vector>& Ef_bound=
        patch().lookupPatchField<volVectorField, scalar>("EField");
    const fvPatchField<scalar>& muef=
        patch().lookupPatchField<volScalarField, scalar>("mobilityelectron");
    const fvPatchField<scalar>& kappae=
        patch().lookupPatchField<volScalarField, scalar>("kappa_electron");
    const volVectorField& gradTe =
        db().objectRegistry::lookupObject<volVectorField>("gradTe");
    const volScalarField& electron =
        db().objectRegistry::lookupObject<volScalarField>("electron");
    const volScalarField& kappae_cell =
        db().objectRegistry::lookupObject<volScalarField>("kappa_electron");
    const fvPatchVectorField& FiT = Fif.boundaryField()[patchi];
    scalarField expterm = FiT&n;
    scalarField pp = FiT&n;
    scalarField qq = FiT&n;
    seec_=0.0;
    forAll(expterm, facei)
    {
        Info <<kappae[facei] <<"Kappa e ===================="<< kappae_cell[facei]<< endl;
    	label faceCelli = patch().faceCells()[facei];
        Info << flux_bound[facei] << flux_cell[faceCelli] << endl;
        int a = pos(Ef[faceCelli]&n[facei]);
        int c = pos(flux_cell[faceCelli]&n[facei]);
        float eteemp = Tef[faceCelli];
        //if (eteemp > 22*11609)
        //    eteemp = 22*11600;
        float vth = sqrt(8.0*1.38e-23*eteemp/9.1e-31/acos(-1.0)) ;
		expterm[facei]=a*(1.38e-23*Tse_*seec_)*(Fif[faceCelli]&n[facei]) + 0*kappae[facei]*(gradTe[faceCelli] & n[facei]);// (2/3)*vth*(1-a)*(1.38e-23*Tse_*seec_)*(Fif[facei]&n[facei])/(muef[facei]*(Ef_bound[facei]&n[facei]))+ 2*(1-a)*(1.38e-23*Tse_*seec_)*(Fif[facei]&n[facei]);
        pp[facei] = c*((5/2)*1.38e-23*(flux_cell[faceCelli]&n[facei]) + vth*1.38e-23*electron[faceCelli])- kappae[facei]*this->patch().deltaCoeffs()[facei];//(1-2*a)*(5/3)*muef[facei]*(Ef[facei]&n[facei]) + (2/3)*vth- (5/3)*kappae[facei]*this->patch().deltaCoeffs()[facei];
        qq[facei] = (5/2)*1.38e-23*(flux_bound[facei]&n[facei]) - kappae[facei]*this->patch().deltaCoeffs()[facei];
        //int a = pos(-Ef[faceCelli]&n[facei]);
        //float vth = sqrt(8.0*1.38e-23*Tef[faceCelli]/9.1e-31/acos(-1.0)) ;
		//expterm[facei]= (2/3)*vth*(1-a)*(1.38e-23*Tse_*seec_)*(Fif[faceCelli]&n[facei])/(muef[facei]*(Ef[faceCelli]&n[facei]))+ 2*(1-a)*(1.38e-23*Tse_*seec_)*(Fif[faceCelli]&n[facei]);
        //pp[facei] =(1-2*a)*(5/3)*muef[facei]*(Ef[faceCelli]&n[facei]) + (2/3)*vth- (5/3)*kappae[facei]*this->patch().deltaCoeffs()[facei];
        //qq[facei] = -(5/3)*muef[facei]*(Ef_bound[faceCelli]&n[facei]) - (5/3)*kappae[facei]*this->patch().deltaCoeffs()[facei];
	}
    //Info << "lets see the value of muef" << muef->patchInternalField() << endl;
    this->valueFraction() = 1-pp/(qq+1e-50); 
    this->refValue() = -expterm/qq/this->valueFraction();
    this->refGrad() = 0;// expterm/pp*this->patch().deltaCoeffs();
    mixedFvPatchField<scalar>::updateCoeffs();
    
}


void Foam::electronTemperatureWilsonBC5::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::electronTemperatureWilsonBC5::operator=
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
        electronTemperatureWilsonBC5
    );
} // End namespace Foam
// ************************************************************************* //
