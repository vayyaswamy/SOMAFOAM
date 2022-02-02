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
    R_(0.0)
    //count_(0.0)
{}


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
    R_(dict.lookupOrDefault<scalar>("R", 0.0))
    //count_(0.0)
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
    amplitude_(ptf.amplitude_, mapper),
    modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    bias_(ptf.bias_),
    R_(ptf.R_)
    //count_(ptf.count_)
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
    R_(tppsf.R_)
    //count_(tppsf.count_)
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
    R_(tppsf.R_)
    //count_(tppsf.count_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::externalCircuitPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volVectorField& Jtot =
        db().objectRegistry::lookupObject<volVectorField>("Jtot");

    label patchi = this->patch().index();
    const fvPatchVectorField& JtotPatch = Jtot.boundaryField()[patchi];

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    scalar patchCurrent = 0.0;

    //Info << "JtotD = " << JtotDPatch << endl;

    forAll (JtotPatch,faceI)
    {
        patchCurrent += JtotPatch[faceI] & mesh.Sf().boundaryField()[patchi][faceI] ;
    }

    Info << "patchCurrent = " << patchCurrent << endl ;

    /*iSqrSum_ += sqr(patchCurrent)*this->db().time().deltaTValue();
    tSum_ += *this->db().time().deltaTValue();

    if (tSum_ >= 1.0/frequency_)
    {

        scalar iRMSold = iRMS_;

        iRMS_ = sqrt(iSqrSum_/tSum_);

        tSum_ = tSum_ - 1.0/frequency_;
        
        iSqrSum_ = sqr(patchCurrent)*tSum_;

        amplitude_ = amplitude_*(iRMS_/iDesired_)
    }*/


    if (modelName_ == "directCurrent")
    {
        Info << "V = " << amplitude_ + patchCurrent*R_;
        operator==(amplitude_ + patchCurrent*R_) ;
    }
    else if (modelName_ == "cosFrequencyModulated")
    {
        operator==(amplitude_*Foam::cos(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_ + patchCurrent*R_);
    }
    else if (modelName_ == "sinFrequencyModulated")
    {
        //Info << "V = " << (amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_ ) << endl;
        //Info << "V = " << (amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_ + patchCurrent*R_) << endl;
        operator==(amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_ + patchCurrent*R_);
    }
    else
    {
        FatalErrorIn
        (
            "externalCircuitPotential::updateCoeffs()"
        )   << " model name inconsitent, model = " << modelName_
            << exit(FatalError);
    }

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
        << R_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, externalCircuitPotential);
}
// ************************************************************************* //
