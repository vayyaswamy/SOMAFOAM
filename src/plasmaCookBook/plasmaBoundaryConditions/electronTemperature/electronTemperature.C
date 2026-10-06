/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electronTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronTemperature::electronTemperature
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    seec_(0),
    Tse_(0),
    Edepend_(false),
    TFN_(0),
    FE_(false),
    beta_(1.0),
    wf_(1.0)
{}


Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
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
    zeroGradientFvPatchScalarField(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Tse_(readScalar(dict.lookup("Tse"))),
    Edepend_(dict.lookupOrDefault<bool>("Edepend", true)),
    TFN_(readScalar(dict.lookup("TFN"))),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::electronTemperature::electronTemperature
(
    const electronTemperature& ptf
)
:
    zeroGradientFvPatchScalarField(ptf),
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
    zeroGradientFvPatchScalarField(ptf, iF),
    seec_(ptf.seec_),
    Tse_(ptf.Tse_),
    Edepend_(ptf.Edepend_),
    TFN_(ptf.TFN_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::electronTemperature::emittedFluxes
(
    scalarField& secondary,
    scalarField& field
) const
{
    const vectorField n(patch().nf());

    const fvPatchField<vector>& Ef =
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<vector>& Fif =
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    const scalarField Enorm(Ef & n);
    const scalarField Fifnorm(Fif & n);

    // Secondary emission by the ions that reach the wall
    secondary = seec_*pos(Fifnorm)*Fifnorm;

    field = 0.0*Enorm;

    if (FE_)
    {
        // Fowler-Nordheim emission where the field points into the wall
        // (as in electronTemperature; 1e-2 converts V/m to V/cm)
        const scalarField c(pos(Enorm));

        const scalarField vofy
        (
            0.95 - sqr(3.79e-4)*beta_*c*mag(Enorm)*1e-2/sqr(wf_)
        );

        field =
            1.54e-6/1.602e-19*sqr(beta_*c*mag(Enorm))/1.1/wf_
           *exp(-6.85e9*pow(wf_, 1.5)*vofy/beta_/(c*mag(Enorm) + SMALL));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::electronTemperature::emittedFlux() const
{
    scalarField secondary;
    scalarField field;

    emittedFluxes(secondary, field);

    return tmp<scalarField>(new scalarField(secondary + field));
}


Foam::tmp<Foam::scalarField>
Foam::electronTemperature::emittedEnergyFlux() const
{
    scalarField secondary;
    scalarField field;

    emittedFluxes(secondary, field);

    return tmp<scalarField>
    (
        new scalarField
        (
            2.0*plasmaConstants::boltzC.value()*(secondary*Tse_ + field*TFN_)
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::electronTemperature::energyFluxPerKelvin
(
    const scalarField& flux
) const
{
    const fvPatchField<scalar>& Nef =
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    const scalarField& Tef = *this;

    // Thermal flux of the plasma electrons to the wall, n vth/4
    const scalarField thermalFlux
    (
        0.25*Nef*sqrt(8.0*1.38e-23*Tef/9.1e-31/acos(-1.0))
    );

    return tmp<scalarField>
    (
        new scalarField
        (
            plasmaConstants::boltzC.value()
           *(
                2.0*min(flux, thermalFlux)
              + 2.5*max(flux - thermalFlux, scalar(0))
            )
        )
    );
}


void Foam::electronTemperature::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("seec") << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tse") << Tse_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend") << Edepend_ << token::END_STATEMENT << nl;
    os.writeKeyword("TFN") << TFN_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_emission") << FE_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_enhancement_factor")
        << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("work_function") << wf_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electronTemperature
    );
}

// ************************************************************************* //
