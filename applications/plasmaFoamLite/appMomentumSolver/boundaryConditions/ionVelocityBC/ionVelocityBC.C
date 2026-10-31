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
#include "ionVelocityBC.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::ionVelocityBC::ionVelocityBC
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    fieldName_(iF.name()),
    mass_(0.0),
    Tion_(0.0)
{
    this->refValue() = vector::zero;
    this->refGrad() = vector::zero;
    this->valueFraction() = 0;
    fvPatchField<vector>::operator=(this->patchInternalField());
}


Foam::ionVelocityBC::ionVelocityBC
(
    const ionVelocityBC& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    mass_(ptf.mass_),
    Tion_(ptf.Tion_)
{
}

Foam::ionVelocityBC::ionVelocityBC
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    fieldName_(iF.name()),
    mass_(readScalar(dict.lookup("mass"))),
    Tion_(readScalar(dict.lookup("Tion")))
{
    this->refValue() = vector::zero;
    this->refGrad() = vector::zero;
    this->valueFraction() = 0;
    fvPatchField<vector>::operator=(this->patchInternalField());
}

Foam::ionVelocityBC::ionVelocityBC
(
    const ionVelocityBC& ptf
)
:
    mixedFvPatchField<vector>(ptf),
    fieldName_(ptf.fieldName_),
    mass_(ptf.mass_),
    Tion_(ptf.Tion_)
{
    fvPatchField<vector>::operator=(this->patchInternalField());
}

Foam::ionVelocityBC::ionVelocityBC
(
    const ionVelocityBC& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    fieldName_(ptf.fieldName_),
    mass_(ptf.mass_),
    Tion_(ptf.Tion_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::ionVelocityBC::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    vectorField n = patch().nf();
    const fvPatchField<vector>& VEf=
        patch().lookupPatchField<volVectorField, vector>(fieldName_);
 
    this->refValue() = 0.25*sqrt(8.0*1.38e-23*Tion_/mass_/acos(-1.0))* n;
    this->valueFraction() = 1-pos(VEf & n) ;
    this->refGrad() = vector::zero;
    mixedFvPatchField<vector>::updateCoeffs();    
}


void Foam::ionVelocityBC::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("mass")
        << mass_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tion")
        << Tion_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::ionVelocityBC::operator=
(
    const fvPatchField<vector>& ptf
)
{
    fvPatchField<vector>::operator=
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
        fvPatchVectorField,
        ionVelocityBC
    );
} // End namespace Foam
// ************************************************************************* //
