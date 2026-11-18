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
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCircuitPotential::externalCircuitPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    amplitude_(p.size(), 0),
    modelName_("directCurrent"),
    frequency_(0.0),
    bias_(0.0),
    R_(0.0),
    C_(1e20),
    rampTime_(0),
    gapCapacitance_(0),
    QC_(0),
    Qe_(0),
    I_(0),
    timeIndex_(-1),
    QeSet_(false),
    lastPhiP_(GREAT)
{}


Foam::externalCircuitPotential::externalCircuitPotential
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
    rampTime_(dict.lookupOrDefault<scalar>("rampTime", 0)),
    gapCapacitance_(dict.lookupOrDefault<scalar>("gapCapacitance", 0)),
    QC_(dict.lookupOrDefault<scalar>("QC", 0)),
    Qe_(dict.lookupOrDefault<scalar>("Qe", 0)),
    I_(dict.lookupOrDefault<scalar>("I", 0)),
    timeIndex_(-1),
    QeSet_(dict.found("Qe")),
    lastPhiP_(GREAT)
{
    if
    (
        modelName_ != "directCurrent"
     && modelName_ != "sinFrequencyModulated"
     && modelName_ != "cosFrequencyModulated"
    )
    {
        FatalIOErrorIn("externalCircuitPotential", dict)
            << "Unknown model " << modelName_ << "; valid models are "
            << "directCurrent, sinFrequencyModulated and "
            << "cosFrequencyModulated" << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}


Foam::externalCircuitPotential::externalCircuitPotential
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
    R_(ptf.R_),
    C_(ptf.C_),
    rampTime_(ptf.rampTime_),
    gapCapacitance_(ptf.gapCapacitance_),
    QC_(ptf.QC_),
    Qe_(ptf.Qe_),
    I_(ptf.I_),
    timeIndex_(ptf.timeIndex_),
    QeSet_(ptf.QeSet_),
    lastPhiP_(ptf.lastPhiP_)
{}


Foam::externalCircuitPotential::externalCircuitPotential
(
    const externalCircuitPotential& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    amplitude_(ptf.amplitude_),
    modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    bias_(ptf.bias_),
    R_(ptf.R_),
    C_(ptf.C_),
    rampTime_(ptf.rampTime_),
    gapCapacitance_(ptf.gapCapacitance_),
    QC_(ptf.QC_),
    Qe_(ptf.Qe_),
    I_(ptf.I_),
    timeIndex_(ptf.timeIndex_),
    QeSet_(ptf.QeSet_),
    lastPhiP_(ptf.lastPhiP_)
{}


Foam::externalCircuitPotential::externalCircuitPotential
(
    const externalCircuitPotential& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    amplitude_(ptf.amplitude_),
    modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    bias_(ptf.bias_),
    R_(ptf.R_),
    C_(ptf.C_),
    rampTime_(ptf.rampTime_),
    gapCapacitance_(ptf.gapCapacitance_),
    QC_(ptf.QC_),
    Qe_(ptf.Qe_),
    I_(ptf.I_),
    timeIndex_(ptf.timeIndex_),
    QeSet_(ptf.QeSet_),
    lastPhiP_(ptf.lastPhiP_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::externalCircuitPotential::electrodeCharge() const
{
    // Q = sum eps dPhi/dn A, n out of the domain: positive on an electrode
    // above the potential of the cells next to it
    const volScalarField& eps =
        db().objectRegistry::lookupObject<volScalarField>("epsilon");

    const scalarField epsP(eps.boundaryField()[patch().index()]);

    return gSum(epsP*snGrad()*patch().magSf());
}


Foam::scalar Foam::externalCircuitPotential::conductionCurrent() const
{
    if (!db().objectRegistry::foundObject<volVectorField>("Jnet"))
    {
        return 0;
    }

    const volVectorField& Jnet =
        db().objectRegistry::lookupObject<volVectorField>("Jnet");

    // From the electrode into the discharge: against the outward normal
    return -gSum(Jnet.boundaryField()[patch().index()] & patch().Sf());
}


Foam::scalar Foam::externalCircuitPotential::gapCapacitance() const
{
    if (gapCapacitance_ > 0)
    {
        return gapCapacitance_;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (mesh.nGeometricD() != 1)
    {
        FatalErrorIn("externalCircuitPotential::gapCapacitance()")
            << "Patch " << patch().name() << ": set gapCapacitance [F], "
            << "the capacitance between the electrode and the other "
            << "electrodes, on a mesh that is not one-dimensional"
            << exit(FatalError);
    }

    // A/sum(dx/eps) along the plasma mesh
    const volScalarField& eps =
        db().objectRegistry::lookupObject<volScalarField>("epsilon");

    const scalar A = gSum(patch().magSf());
    const scalarField& V = mesh.V();

    scalar sumDxEps = 0;
    forAll(V, cellI)
    {
        sumDxEps += V[cellI]/A/eps[cellI];
    }
    reduce(sumDxEps, sumOp<scalar>());

    return A/sumDxEps;
}


Foam::scalar Foam::externalCircuitPotential::sourceVoltage
(
    const scalar t
) const
{
    const scalar a = gAverage(amplitude_);

    if (modelName_ == "directCurrent")
    {
        const scalar ramp =
            (rampTime_ > 0) ? 1 - Foam::exp(-t/rampTime_) : 1;

        return a*ramp + bias_;
    }
    else if (modelName_ == "sinFrequencyModulated")
    {
        return a*Foam::sin(2*mathematicalConstant::pi*frequency_*t) + bias_;
    }

    return a*Foam::cos(2*mathematicalConstant::pi*frequency_*t) + bias_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCircuitPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const Time& runTime = db().time();
    const scalar t = runTime.value();

    if (modelName_ == "cosFrequencyModulated")
    {
        operator==(sourceVoltage(t));
        fixedValueFvPatchScalarField::updateCoeffs();
        return;
    }

    // Solve the circuit once per Poisson solve: if the field next to the
    // electrode has not changed since the last evaluation in this step, the
    // charge would mix the new patch potential with the old field
    const scalar PhiP = gSum(patchInternalField());

    if (runTime.timeIndex() == timeIndex_ && PhiP == lastPhiP_)
    {
        fixedValueFvPatchScalarField::updateCoeffs();
        return;
    }

    lastPhiP_ = PhiP;

    // Present iterate: potential and charge of the electrode
    const scalar Vstar = gAverage(*this);
    const scalar Qstar = electrodeCharge();

    // First evaluation in a new time step: the field holds the solution of
    // the last step, whose charge and current close that step
    if (runTime.timeIndex() != timeIndex_)
    {
        if (QeSet_ && timeIndex_ >= 0)
        {
            const scalar dt0 = runTime.deltaT0Value();
            I_ = (Qstar - Qe_)/dt0 + conductionCurrent();
            QC_ += I_*dt0;
        }

        Qe_ = Qstar;
        QeSet_ = true;
        timeIndex_ = runTime.timeIndex();
    }

    const scalar dt = runTime.deltaTValue();
    const scalar Cg = gapCapacitance();
    const scalar Reff = R_ + dt/C_;
    const scalar Ic = conductionCurrent();

    const scalar V =
        (
            sourceVoltage(t) - QC_/C_
          - Reff*((Qstar - Cg*Vstar - Qe_)/dt + Ic)
        )
       /(1 + Reff*Cg/dt);

    // Current of this iterate, for output
    I_ = (Qstar + Cg*(V - Vstar) - Qe_)/dt + Ic;

    operator==(V);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::externalCircuitPotential::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("model") << modelName_ << token::END_STATEMENT << nl;
    amplitude_.writeEntry("amplitude", os);
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("bias") << bias_ << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
    os.writeKeyword("rampTime") << rampTime_ << token::END_STATEMENT << nl;
    if (gapCapacitance_ > 0)
    {
        os.writeKeyword("gapCapacitance")
            << gapCapacitance_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("QC") << QC_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qe") << Qe_ << token::END_STATEMENT << nl;
    os.writeKeyword("I") << I_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, externalCircuitPotential);
}

// ************************************************************************* //
