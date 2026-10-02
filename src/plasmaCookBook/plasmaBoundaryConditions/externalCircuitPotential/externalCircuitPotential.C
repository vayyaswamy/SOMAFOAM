/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "externalCircuitPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mathematicalConstants.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::externalCircuitPotential::
externalCircuitPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    amplitude_(p.size(), 0.0),
    modelName_("directCurrent"),
    frequency_(0.0),
    bias_(0.0),
    R_(0.0),
    C_(1e20),
    Phi_(p.size(), 0),
    I_(0)
{
}


Foam::externalCircuitPotential::
externalCircuitPotential
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
    bias_(dict.lookupOrDefault<scalar>("bias", 0.0)),
    R_(dict.lookupOrDefault<scalar>("R", 0.0)),
    C_(dict.lookupOrDefault<scalar>("C", 1e20)),
    Phi_("Phi", dict, p.size()),
    I_(dict.lookupOrDefault<scalar>("I", 0))
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

Foam::externalCircuitPotential::
externalCircuitPotential
(
    const externalCircuitPotential& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    bias_(ptf.bias_),
    R_(ptf.R_),
    C_(ptf.C_),
    Phi_(ptf.Phi_),
    I_(ptf.I_)
{}


Foam::externalCircuitPotential::
externalCircuitPotential
(
    const externalCircuitPotential& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_),
    R_(tppsf.R_),
    C_(tppsf.C_),
    Phi_(tppsf.Phi_),
    I_(tppsf.I_)
{}


Foam::externalCircuitPotential::
externalCircuitPotential
(
    const externalCircuitPotential& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_),
    R_(tppsf.R_),
    C_(tppsf.C_),
    Phi_(tppsf.Phi_),
    I_(tppsf.I_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::externalCircuitPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volScalarField& QC =
        db().objectRegistry::lookupObject<volScalarField>("QC");

    const volVectorField& Jtotpatch =
        db().objectRegistry::lookupObject<volVectorField>("Jtot");

    label patchi = this->patch().index();
    const fvPatchScalarField& QClocal = QC.boundaryField()[patchi];

    const fvMesh& mesh = patch().boundaryMesh().mesh();


    scalar patchCurrent = gSum(Jtotpatch & mesh.Sf().boundaryField()[patchi]);

    Info << "patchCurrent = " << patchCurrent << endl;


    if (modelName_ == "directCurrent")
    {
        Info << "V = " << amplitude_*(1-Foam::exp(-this->db().time().value()/C_)) + neg(patchCurrent)*patchCurrent*R_;
        operator==(amplitude_*(1-Foam::exp(-this->db().time().value()/C_)) + neg(patchCurrent)*patchCurrent*R_) ;
    }
    else if (modelName_ == "cosFrequencyModulated")
    {
        operator==(amplitude_*Foam::cos(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_);
    }
    else if (modelName_ == "sinFrequencyModulated")
    {
        Phi_ = amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_ + patchCurrent*R_ + QClocal/C_;
        Info << "QC = " << QClocal[0] << endl;
        Info << "Applied = " << amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) << endl;
        operator==(Phi_);
    }
    else
    {
        FatalErrorIn
        (
            "externalCircuitPotential::updateCoeffs()"
        )   << " model name inconsitent, model = " << modelName_
            << exit(FatalError);
    }
    Info << "Phi_ = " << Phi_ << endl;
    Info << "I_ = " << I_ << endl;
    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::externalCircuitPotential::
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
    os.writeKeyword("R")
        << R_  << token::END_STATEMENT << nl;
    os.writeKeyword("C")
        << C_ << token::END_STATEMENT << nl;
    Phi_.writeEntry("Phi", os);
    os.writeKeyword("I")
        << I_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, externalCircuitPotential);
}
// ************************************************************************* //
