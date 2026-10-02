/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electronMomentumTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronMomentumTemperature::electronMomentumTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(0),
    Tse_(0),
    TFN_(0),
    Edepend_(true),
    driftEnergyFactor_(2.5),
    FE_(false),
    beta_(1.0),
    wf_(1.0)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
}


Foam::electronMomentumTemperature::electronMomentumTemperature
(
    const electronMomentumTemperature& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    TFN_(ptf.TFN_),
    Edepend_(ptf.Edepend_),
    driftEnergyFactor_(ptf.driftEnergyFactor_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


Foam::electronMomentumTemperature::electronMomentumTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Tse_(readScalar(dict.lookup("Tse"))),
    TFN_(readScalar(dict.lookup("TFN"))),
    Edepend_(dict.lookupOrDefault<Switch>("Edepend", true)),
    driftEnergyFactor_(dict.lookupOrDefault<scalar>("driftEnergyFactor", 2.5)),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;

    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::electronMomentumTemperature::electronMomentumTemperature
(
    const electronMomentumTemperature& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    TFN_(ptf.TFN_),
    Edepend_(ptf.Edepend_),
    driftEnergyFactor_(ptf.driftEnergyFactor_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


Foam::electronMomentumTemperature::electronMomentumTemperature
(
    const electronMomentumTemperature& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    TFN_(ptf.TFN_),
    Edepend_(ptf.Edepend_),
    driftEnergyFactor_(ptf.driftEnergyFactor_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electronMomentumTemperature::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar kB = 1.38e-23;

    const vectorField n = patch().nf();

    const fvPatchField<scalar>& Nef =
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    const fvPatchField<vector>& Uef =
        patch().lookupPatchField<volVectorField, vector>("U_electron");

    const fvPatchField<scalar>& kappaef =
        patch().lookupPatchField<volScalarField, scalar>("kappa_electron");

    const fvPatchField<vector>& Ef =
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<vector>& Fif =
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    // Electrons injected by the wall (same models as electronTemperature)
    const scalarField Fifnorm = Fif & n;

    const scalarField Gamma_se = seec_*pos(Fifnorm)*Fifnorm;

    scalarField Gamma_FE(this->size(), 0.0);

    if (FE_)
    {
        const scalarField Enorm = Ef & n;

        // only if the E-field points into the wall
        const scalarField c = pos(Enorm);

        // 1E-2 converts V/m to V/cm
        const scalarField vofy =
            0.95 - sqr(3.79E-4)*beta_*c*mag(Enorm)*1E-2/sqr(wf_);

        Gamma_FE =
            1.54E-6/1.602e-19*sqr(beta_*c*mag(Enorm))/1.1/wf_
           *exp(-6.85E9*pow(wf_, 1.5)*vofy/beta_/(c*mag(Enorm) + SMALL));
    }

    const scalarField Gamma_inj = Gamma_se + Gamma_FE;

    // Net electron wall flux as imposed by the electron velocity condition,
    // and the part of it leaving through the wall
    const scalarField Gamma = Nef*(Uef & n);

    const scalarField Gamma_out = max(Gamma + Gamma_inj, scalar(0));

    // Energy per outgoing electron: thermal loss (2 k Te) where the field
    // repels electrons, drift loss where it drives them into the wall
    scalarField a(this->size(), 1.0);

    if (Edepend_)
    {
        a = pos(Ef & n);
    }

    const scalarField eps = a*2.0 + (1.0 - a)*driftEnergyFactor_;

    // kappa dTe/dn = C1 Te + C3
    const scalarField C1 = kB*((2.5 - eps)*Gamma_out - 2.5*Gamma_inj);

    const scalarField C2 = kappaef;

    const scalarField C3 = kB*(Tse_*Gamma_se + TFN_*Gamma_FE);

    this->refValue() = 0.0;

    this->valueFraction() = C1/(C1 - C2*this->patch().deltaCoeffs());

    this->refGrad() = C3/(C2 + SMALL);

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::electronMomentumTemperature::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse")
        << Tse_ << token::END_STATEMENT << nl;
    os.writeKeyword("TFN")
        << TFN_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend")
        << Edepend_ << token::END_STATEMENT << nl;
    os.writeKeyword("driftEnergyFactor")
        << driftEnergyFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_emission")
        << FE_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_enhancement_factor")
        << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("work_function")
        << wf_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::electronMomentumTemperature::operator=
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
        electronMomentumTemperature
    );
}

// ************************************************************************* //
