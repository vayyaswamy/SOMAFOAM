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

#include "electronTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"
#include "surfaceFields.H"
#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronTemperature::electronTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0),
    Tse_(0),
    Edepend_(false),
    TFN_(0),
    FE_(false),
    beta_(1.0),
    wf_(1.0)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;

}


Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_),
    TFN_(ptf.TFN_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


Foam::electronTemperature::electronTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Tse_(readScalar(dict.lookup("Tse"))),
    Edepend_(readBool(dict.lookup("Edepend"))),
    TFN_(readScalar(dict.lookup("TFN"))),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
{
    this->refValue() = 0;

    this->refGrad() = 0;
    this->valueFraction() = 0;

    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_),
    TFN_(ptf.TFN_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_),
    TFN_(ptf.TFN_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::electronTemperature::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const vectorField n(patch().nf());

    const fvPatchField<scalar>& kappaef =
        patch().lookupPatchField<volScalarField, scalar>("kappa_electron");

    const fvPatchField<vector>& Ef =
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<vector>& Fif =
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    const scalarField Enorm(Ef & n);
    const scalarField Fifnorm(Fif & n);

    // Emitted electrons: secondaries from the ions that reach the wall and
    // field emission where the field points into the wall
    const scalarField Gamma_se(seec_*pos(Fifnorm)*Fifnorm);

    scalarField Gamma_FE(patch().size(), 0.0);

    if (FE_)
    {
        const scalarField c(pos(Enorm));

        // 1e-2 converts V/m to V/cm
        const scalarField vofy
        (
            0.95 - sqr(3.79e-4)*beta_*c*mag(Enorm)*1e-2/sqr(wf_)
        );

        Gamma_FE =
            1.54e-6/1.602e-19*sqr(beta_*c*mag(Enorm))/1.1/wf_
           *exp(-6.85e9*pow(wf_, 1.5)*vofy/beta_/(c*mag(Enorm) + SMALL));
    }

    // Energy flux to the wall (Lymberopoulos and Economou):
    //     q.n = 5/2 k Te Gamma_out - 5/2 k (Gamma_se Tse + Gamma_FE TFN),
    // with Gamma_out the flux of plasma electrons to the wall. The energy
    // equation carries 5/2 k Te Gamma_net - kappa dTe/dn through the wall
    // face, and Gamma_net = Gamma_out - Gamma_se - Gamma_FE, so that
    //     kappa dTe/dn = -hEmitted (Te - Temitted),
    // with hEmitted = 5/2 k (Gamma_se + Gamma_FE) and Temitted the mean
    // temperature of the emitted electrons. hEmitted is never negative:
    // the weight of the condition lies between 0 and 1 for any wall cell,
    // and without emission it is a zero gradient
    const scalarField emitted(Gamma_se + Gamma_FE);

    const scalarField hEmitted(2.5*plasmaConstants::boltzC.value()*emitted);

    this->refValue() =
        (Gamma_se*Tse_ + Gamma_FE*TFN_ + VSMALL*Tse_)/(emitted + VSMALL);

    this->refGrad() = 0.0;

    this->valueFraction() =
        hEmitted/(hEmitted + kappaef*this->patch().deltaCoeffs() + VSMALL);

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::electronTemperature::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend")
        << Edepend_ << token::END_STATEMENT << nl;
    os.writeKeyword("TFN")
        << TFN_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_emission")
        << FE_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_enhancement_factor")
        << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("work_function")
        << wf_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::electronTemperature::operator=
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
        electronTemperature
    );
} // End namespace Foam

// ************************************************************************* //
