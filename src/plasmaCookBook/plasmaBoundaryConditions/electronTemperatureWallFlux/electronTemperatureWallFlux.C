/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electronTemperatureWallFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronTemperatureWallFlux::electronTemperatureWallFlux
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


Foam::electronTemperatureWallFlux::electronTemperatureWallFlux
(
    const electronTemperatureWallFlux& ptf,
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


Foam::electronTemperatureWallFlux::electronTemperatureWallFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Tse_(readScalar(dict.lookup("Tse"))),
    Edepend_(readBool(dict.lookup("Edepend"))),
    TFN_(readScalar(dict.lookup("TFN"))),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
{
    fvPatchField<scalar>::operator=(this->patchInternalField());
}


Foam::electronTemperatureWallFlux::electronTemperatureWallFlux
(
    const electronTemperatureWallFlux& ptf
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


Foam::electronTemperatureWallFlux::electronTemperatureWallFlux
(
    const electronTemperatureWallFlux& ptf,
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

void Foam::electronTemperatureWallFlux::emittedFluxes
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
Foam::electronTemperatureWallFlux::emittedFlux() const
{
    scalarField secondary;
    scalarField field;

    emittedFluxes(secondary, field);

    return tmp<scalarField>(new scalarField(secondary + field));
}


Foam::tmp<Foam::scalarField>
Foam::electronTemperatureWallFlux::emittedEnergyFlux() const
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
Foam::electronTemperatureWallFlux::energyPerElectron() const
{
    tmp<scalarField> tw
    (
        new scalarField(patch().size(), 2.0*plasmaConstants::boltzC.value())
    );

    if (Edepend_)
    {
        const fvPatchField<vector>& Ef =
            patch().lookupPatchField<volVectorField, vector>("E");

        // 2 k where the electrons leave by their thermal motion (field
        // pointing into the wall), 5/2 k where they drift to the wall
        const scalarField thermal(pos(Ef & patch().nf()));

        tw() = plasmaConstants::boltzC.value()*(2.5 - 0.5*thermal);
    }

    return tw;
}


void Foam::electronTemperatureWallFlux::write(Ostream& os) const
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
        electronTemperatureWallFlux
    );
}

// ************************************************************************* //
