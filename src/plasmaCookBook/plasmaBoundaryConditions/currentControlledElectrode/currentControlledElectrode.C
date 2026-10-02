/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "currentControlledElectrode.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mathematicalConstants.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::currentControlledElectrode::
currentControlledElectrode
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    amplitude_(p.size(), 0.0),
    modelName_("directCurrent"),
    frequency_(0.0),
    bias_(0.0)
{}


Foam::currentControlledElectrode::
currentControlledElectrode
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    amplitude_("amplitude", dict, p.size()),
    modelName_(dict.lookupOrDefault<word>("model", "directCurrent")),
    frequency_(dict.lookupOrDefault<scalar>("frequency", 0.0)),
    bias_(dict.lookupOrDefault<scalar>("bias", 0.0))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}

Foam::currentControlledElectrode::
currentControlledElectrode
(
    const currentControlledElectrode& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_, mapper),
    modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    bias_(ptf.bias_)
{}


Foam::currentControlledElectrode::
currentControlledElectrode
(
    const currentControlledElectrode& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_)
{}


Foam::currentControlledElectrode::
currentControlledElectrode
(
    const currentControlledElectrode& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::currentControlledElectrode::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    amplitude_.autoMap(m);
    modelName_.autoMap(m);
    frequency_.autoMap(m);
    bias_.autoMap(m);
}


void Foam::currentControlledElectrode::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const currentControlledElectrode& tiptf =
        refCast<const currentControlledElectrode>(ptf);

    amplitude_.rmap(tiptf.amplitude_, addr);
    modelName_.rmap(tiptf.modelName_, addr);
    frequency_.rmap(tiptf.frequency_, addr);
    bias_.rmap(tiptf.bias_, addr);
}

void Foam::currentControlledElectrode::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatch& patch = this->patch();

    if (modelName_ == "directCurrent")
    {
        operator==(amplitude_);
    }
    else if (modelName_ == "cosFrequencyModulated")
    {
        scalar desired_current_density = (amplitude_*Foam::cos(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_);
        const fvPatchField<vector>& Ef = patch().lookupPatchField<volVectorField, vector>("E");
        const fvPatchField<vector>& Jnet = patch().lookupPatchField<volVectorField, vector>("Jnet");
    }
    else if (modelName_ == "sinFrequencyModulated")
    {
        operator==(amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_);
    }
    else
    {
        FatalErrorIn
        (
            "plasmaPotential::updateCoeffs()"
        )   << " model name inconsitent, model = " << modelName_
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::plasmaPotential::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("model")
        << modelName_ << token::END_STATEMENT << nl;
    amplitude_.writeEntry("amplitude", os);
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("bias")
        << bias_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, plasmaPotential);
}
// ************************************************************************* //
