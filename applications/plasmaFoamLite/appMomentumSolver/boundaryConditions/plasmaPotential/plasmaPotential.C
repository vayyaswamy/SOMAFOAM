/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors
(This code was written by Abhishek Kumar Verma a former PhD student at UC Merced)
License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mathematicalConstants.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::plasmaPotential::
plasmaPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    amplitude_(p.size(), 0.0),
    modelName_("directCurrent"),
    frequency_(0.0),
    bias_(0.0)
{}


Foam::plasmaPotential::
plasmaPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
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

Foam::plasmaPotential::
plasmaPotential
(
    const plasmaPotential& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_, mapper),
    modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    bias_(ptf.bias_)
{}


Foam::plasmaPotential::
plasmaPotential
(
    const plasmaPotential& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_)
{}


Foam::plasmaPotential::
plasmaPotential
(
    const plasmaPotential& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaPotential::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    amplitude_.autoMap(m);
}


void Foam::plasmaPotential::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const plasmaPotential& tiptf =
        refCast<const plasmaPotential>(ptf);

    amplitude_.rmap(tiptf.amplitude_, addr);
}
void Foam::plasmaPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (modelName_ == "directCurrent")
    {
        operator==(amplitude_);
    }
    else if (modelName_ == "sinFrequencyModulated")
    {
        operator==(amplitude_*Foam::cos(2*22/7*frequency_*this->db().time().value()) + bias_);
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
