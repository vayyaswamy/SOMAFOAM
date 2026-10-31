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
#include "zeroFluxElectronDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroFluxElectronDensity::zeroFluxElectronDensity
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

Foam::zeroFluxElectronDensity::zeroFluxElectronDensity
(
    const zeroFluxElectronDensity& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_)
{
}


Foam::zeroFluxElectronDensity::zeroFluxElectronDensity
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



Foam::zeroFluxElectronDensity::zeroFluxElectronDensity
(
    const zeroFluxElectronDensity& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_)
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}

Foam::zeroFluxElectronDensity::zeroFluxElectronDensity
(
    const zeroFluxElectronDensity& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::zeroFluxElectronDensity::updateCoeffs()
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
    forAll(temp1, facei)
    {
        temp1[facei]=(velocityCell[facei]&n[facei])/((vel_bound[facei]&n[facei])+1e-50);

	}
    //Info << "lets see the value of muef" << muef->patchInternalField() << endl;
     this->refValue() = 0;
    this->valueFraction() = 1-temp1;
    this->refGrad() = 0;
    mixedFvPatchField<scalar>::updateCoeffs();
}

void Foam::zeroFluxElectronDensity::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void Foam::zeroFluxElectronDensity::operator=
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
        zeroFluxElectronDensity
    );
} // End namespace Foam
// ************************************************************************* //
