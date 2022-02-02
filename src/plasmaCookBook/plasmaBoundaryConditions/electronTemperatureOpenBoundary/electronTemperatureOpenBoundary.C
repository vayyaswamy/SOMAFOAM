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

#include "electronTemperatureOpenBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronTemperatureOpenBoundary::electronTemperatureOpenBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    //Info << "Constructor 1" << endl;
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::electronTemperatureOpenBoundary::electronTemperatureOpenBoundary
(
    const electronTemperatureOpenBoundary& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::electronTemperatureOpenBoundary::electronTemperatureOpenBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    //Info << "Constructor 2 " << endl;
    //Info << "fieldName = " << iF.name() << endl;
    this->refValue() = 0.0;

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
    fvPatchField<scalar>::operator=(this->patchInternalField());
}



Foam::electronTemperatureOpenBoundary::electronTemperatureOpenBoundary
(
    const electronTemperatureOpenBoundary& ptf
)
:
    mixedFvPatchField<scalar>(ptf)
{}



Foam::electronTemperatureOpenBoundary::electronTemperatureOpenBoundary
(
    const electronTemperatureOpenBoundary& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::electronTemperatureOpenBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Info << "Inside updateCoeffs " << endl;
    vectorField n = patch().nf();

    

    const fvPatchField<scalar>& Nef=
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    const fvPatchField<scalar>& kappaef=
        patch().lookupPatchField<volScalarField, scalar>("kappa_electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te");

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const scalarField C1 = -chargeSign_*muf*(Ef&n);

    const scalarField C2 = Df;

    scalarField Enorm = Ef&n ;

    scalarField a = pos(Enorm*chargeSign_);

    //Info << "Ef = " << Ef << endl;

    //Info << "n = " << n << endl;

    //Info << "Enorm = " << Enorm << endl;

    //Info << "a = " << a << endl;
 
    this->refValue() = 0.0;

    this->valueFraction() = (1-a)*C1/(C1+C2*this->patch().deltaCoeffs());

    //this->valueFraction() = 0.0;

    this->refGrad() = 0;

    mixedFvPatchField<scalar>::updateCoeffs();

    //Info << "Nef = " << Nef << endl;
}


void Foam::driftDiffusionOpenBoundary::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("chargeSign")
        << chargeSign_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::driftDiffusionOpenBoundary::operator=
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
        driftDiffusionOpenBoundary
    );
} // End namespace Foam

// ************************************************************************* //
