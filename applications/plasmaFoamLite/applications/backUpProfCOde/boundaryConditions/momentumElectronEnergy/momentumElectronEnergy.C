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
#include "momentumElectronEnergy.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::momentumElectronEnergy::momentumElectronEnergy
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


Foam::momentumElectronEnergy::momentumElectronEnergy
(
    const momentumElectronEnergy& ptf,
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

Foam::momentumElectronEnergy::momentumElectronEnergy
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

Foam::momentumElectronEnergy::momentumElectronEnergy
(
    const momentumElectronEnergy& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::momentumElectronEnergy::momentumElectronEnergy
(
    const momentumElectronEnergy& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::momentumElectronEnergy::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
 vectorField n = patch().nf();//normal vector
    label patchi = this->patch().index();
    const volVectorField& velocityCell =
        db().objectRegistry::lookupObject<volVectorField>("velocityelectron");
    const volScalarField& Tef =
        db().objectRegistry::lookupObject<volScalarField>("eTemp");
    const volVectorField& Fif =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux"); 

    const fvPatchField<vector>& vel_bound=
        patch().lookupPatchField<volVectorField, scalar>("velocityelectron");
    const fvPatchVectorField& FiT = Fif.boundaryField()[patchi];
    scalarField temp1 = FiT&n;
    scalarField temp2 = FiT&n;
    //Def= 0*Def;
    forAll(temp1, facei)
    {
    	label faceCelli = patch().faceCells()[facei];
        int a = pos(velocityCell[facei]&n[facei]);
        float vth = sqrt(8.0*1.38e-23*Tef[facei]/9.1e-31/acos(-1.0));
        temp1[facei]=(0.25*vth  + a*(5/3)*(velocityCell[facei]&n[facei]));
        temp2[facei]=(5/2)*1.38e-23*Tse_*a*seec_*(Fif[facei]&n[facei])*this->patch().deltaCoeffs()[facei];
	}
    this->refValue() = 0;
    this->valueFraction() = 1-temp1/(((5/3)*vel_bound&n)+1e-38);
    this->refGrad() = temp2/(temp1+1e-38);
    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::momentumElectronEnergy::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::momentumElectronEnergy::operator=
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
        momentumElectronEnergy
    );
} // End namespace Foam
// ************************************************************************* //
