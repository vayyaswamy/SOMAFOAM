/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "rfPowerElectrode.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mathematicalConstants.H"
#include "foamTime.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rfPowerElectrode::rfPowerElectrode
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    frequency_(0),
    power_(0),
    updateCycles_(5),
    relaxation_(0.3),
    maxChange_(0.1),
    blockingCapacitor_("ideal"),
    capacitance_(1e-12),
    area_(-1),
    amplitude_(0),
    bias_(0),
    charge_(0),
    lastVoltage_(0),
    cycle_(-1),
    cycleTime_(0),
    cycleCharge_(0),
    cycleEnergy_(0),
    updatePeriods_(0),
    updateTime_(0),
    updateEnergy_(0),
    curTimeIndex_(-1)
{}


Foam::rfPowerElectrode::rfPowerElectrode
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    frequency_(readScalar(dict.lookup("frequency"))),
    power_(readScalar(dict.lookup("power"))),
    updateCycles_(dict.lookupOrDefault<label>("updateCycles", 5)),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 0.3)),
    maxChange_(dict.lookupOrDefault<scalar>("maxChange", 0.1)),
    blockingCapacitor_(dict.lookupOrDefault<word>("blockingCapacitor", "ideal")),
    capacitance_(dict.lookupOrDefault<scalar>("capacitance", 1e-12)),
    area_(dict.lookupOrDefault<scalar>("area", -1)),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    bias_(dict.lookupOrDefault<scalar>("bias", 0)),
    charge_(dict.lookupOrDefault<scalar>("charge", 0)),
    lastVoltage_(dict.lookupOrDefault<scalar>("lastVoltage", 0)),
    cycle_(dict.lookupOrDefault<label>("cycle", -1)),
    cycleTime_(dict.lookupOrDefault<scalar>("cycleTime", 0)),
    cycleCharge_(dict.lookupOrDefault<scalar>("cycleCharge", 0)),
    cycleEnergy_(dict.lookupOrDefault<scalar>("cycleEnergy", 0)),
    updatePeriods_(dict.lookupOrDefault<label>("updatePeriods", 0)),
    updateTime_(dict.lookupOrDefault<scalar>("updateTime", 0)),
    updateEnergy_(dict.lookupOrDefault<scalar>("updateEnergy", 0)),
    curTimeIndex_(-1)
{
    if
    (
        blockingCapacitor_ != "ideal"
     && blockingCapacitor_ != "physical"
     && blockingCapacitor_ != "none"
    )
    {
        FatalIOErrorIn
        (
            "rfPowerElectrode::rfPowerElectrode(...)",
            dict
        )   << "Unknown blockingCapacitor " << blockingCapacitor_
            << "; valid entries are ideal, physical and none"
            << exit(FatalIOError);
    }

    if
    (
        frequency_ <= 0 || power_ <= 0 || updateCycles_ < 1
     || relaxation_ <= 0 || relaxation_ > 1 || maxChange_ <= 0
     || capacitance_ <= 0
    )
    {
        FatalIOErrorIn
        (
            "rfPowerElectrode::rfPowerElectrode(...)",
            dict
        )   << "Need frequency > 0, power > 0, updateCycles >= 1,"
            << " 0 < relaxation <= 1, maxChange > 0 and capacitance > 0"
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(bias_);
    }
}


Foam::rfPowerElectrode::rfPowerElectrode
(
    const rfPowerElectrode& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    frequency_(ptf.frequency_),
    power_(ptf.power_),
    updateCycles_(ptf.updateCycles_),
    relaxation_(ptf.relaxation_),
    maxChange_(ptf.maxChange_),
    blockingCapacitor_(ptf.blockingCapacitor_),
    capacitance_(ptf.capacitance_),
    area_(ptf.area_),
    amplitude_(ptf.amplitude_),
    bias_(ptf.bias_),
    charge_(ptf.charge_),
    lastVoltage_(ptf.lastVoltage_),
    cycle_(ptf.cycle_),
    cycleTime_(ptf.cycleTime_),
    cycleCharge_(ptf.cycleCharge_),
    cycleEnergy_(ptf.cycleEnergy_),
    updatePeriods_(ptf.updatePeriods_),
    updateTime_(ptf.updateTime_),
    updateEnergy_(ptf.updateEnergy_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


Foam::rfPowerElectrode::rfPowerElectrode(const rfPowerElectrode& ptf)
:
    fixedValueFvPatchScalarField(ptf),
    frequency_(ptf.frequency_),
    power_(ptf.power_),
    updateCycles_(ptf.updateCycles_),
    relaxation_(ptf.relaxation_),
    maxChange_(ptf.maxChange_),
    blockingCapacitor_(ptf.blockingCapacitor_),
    capacitance_(ptf.capacitance_),
    area_(ptf.area_),
    amplitude_(ptf.amplitude_),
    bias_(ptf.bias_),
    charge_(ptf.charge_),
    lastVoltage_(ptf.lastVoltage_),
    cycle_(ptf.cycle_),
    cycleTime_(ptf.cycleTime_),
    cycleCharge_(ptf.cycleCharge_),
    cycleEnergy_(ptf.cycleEnergy_),
    updatePeriods_(ptf.updatePeriods_),
    updateTime_(ptf.updateTime_),
    updateEnergy_(ptf.updateEnergy_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


Foam::rfPowerElectrode::rfPowerElectrode
(
    const rfPowerElectrode& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    frequency_(ptf.frequency_),
    power_(ptf.power_),
    updateCycles_(ptf.updateCycles_),
    relaxation_(ptf.relaxation_),
    maxChange_(ptf.maxChange_),
    blockingCapacitor_(ptf.blockingCapacitor_),
    capacitance_(ptf.capacitance_),
    area_(ptf.area_),
    amplitude_(ptf.amplitude_),
    bias_(ptf.bias_),
    charge_(ptf.charge_),
    lastVoltage_(ptf.lastVoltage_),
    cycle_(ptf.cycle_),
    cycleTime_(ptf.cycleTime_),
    cycleCharge_(ptf.cycleCharge_),
    cycleEnergy_(ptf.cycleEnergy_),
    updatePeriods_(ptf.updatePeriods_),
    updateTime_(ptf.updateTime_),
    updateEnergy_(ptf.updateEnergy_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rfPowerElectrode::advance()
{
    const Time& runTime = this->db().time();

    // The current and the voltage of the time step that has just been
    // completed: Jtot was evaluated at its end, with the voltage lastVoltage
    const scalar tLast = runTime.value() - runTime.deltaT().value();
    const scalar dtLast = runTime.deltaT0().value();

    if (cycle_ < 0)
    {
        // First call: nothing to add yet
        cycle_ = label(floor(tLast*frequency_ + 1e-9));
        return;
    }

    const volVectorField& Jtot =
        db().objectRegistry::lookupObject<volVectorField>("Jtot");

    const vectorField& Sf = patch().Sf();

    scalar scale = 1;

    if (area_ > 0)
    {
        scale = area_/gSum(patch().magSf());
    }

    // Total current from the electrode into the plasma
    const scalar current =
        -scale*gSum(Jtot.boundaryField()[patch().index()] & Sf);

    cycleTime_ += dtLast;
    cycleCharge_ += current*dtLast;
    cycleEnergy_ += lastVoltage_*current*dtLast;

    if (blockingCapacitor_ == "physical")
    {
        charge_ += current*dtLast;
    }

    // Has the time step that was added completed a period?
    const label cycleNow = label(floor(tLast*frequency_ + 1e-9));

    if (cycleNow <= cycle_)
    {
        return;
    }

    cycle_ = cycleNow;

    const scalar meanCurrent = cycleCharge_/cycleTime_;
    const scalar meanPower = cycleEnergy_/cycleTime_;

    updatePeriods_++;
    updateTime_ += cycleTime_;
    updateEnergy_ += cycleEnergy_;

    if (blockingCapacitor_ == "ideal")
    {
        // The mean current of the period charges the capacitor
        bias_ -= meanCurrent/frequency_/capacitance_;
    }

    if (updatePeriods_ >= updateCycles_)
    {
        const scalar measured = updateEnergy_/updateTime_;

        scalar factor = 1 + maxChange_;

        if (measured > 0)
        {
            factor = pow(power_/measured, 0.5*relaxation_);
            factor = max(min(factor, 1 + maxChange_), 1/(1 + maxChange_));
        }

        amplitude_ *= factor;

        updatePeriods_ = 0;
        updateTime_ = 0;
        updateEnergy_ = 0;
    }

    Info<< "rfPowerElectrode " << patch().name()
        << ": t " << tLast
        << " A " << amplitude_
        << " Vdc "
        << (blockingCapacitor_ == "physical" ? bias_ - charge_/capacitance_ : bias_)
        << " P " << meanPower
        << " I " << meanCurrent << endl;

    cycleTime_ = 0;
    cycleCharge_ = 0;
    cycleEnergy_ = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rfPowerElectrode::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const Time& runTime = this->db().time();

    if (curTimeIndex_ != runTime.timeIndex())
    {
        advance();
        curTimeIndex_ = runTime.timeIndex();
    }

    scalar Vdc = bias_;

    if (blockingCapacitor_ == "physical")
    {
        Vdc -= charge_/capacitance_;
    }

    lastVoltage_ =
        amplitude_
       *Foam::sin(2*mathematicalConstant::pi*frequency_*runTime.value())
      + Vdc;

    operator==(lastVoltage_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::rfPowerElectrode::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("power") << power_ << token::END_STATEMENT << nl;
    os.writeKeyword("updateCycles")
        << updateCycles_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxation") << relaxation_ << token::END_STATEMENT << nl;
    os.writeKeyword("maxChange") << maxChange_ << token::END_STATEMENT << nl;
    os.writeKeyword("blockingCapacitor")
        << blockingCapacitor_ << token::END_STATEMENT << nl;
    os.writeKeyword("capacitance")
        << capacitance_ << token::END_STATEMENT << nl;

    if (area_ > 0)
    {
        os.writeKeyword("area") << area_ << token::END_STATEMENT << nl;
    }

    // State, for restarts
    os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("bias") << bias_ << token::END_STATEMENT << nl;
    os.writeKeyword("charge") << charge_ << token::END_STATEMENT << nl;
    os.writeKeyword("lastVoltage")
        << lastVoltage_ << token::END_STATEMENT << nl;
    os.writeKeyword("cycle") << cycle_ << token::END_STATEMENT << nl;
    os.writeKeyword("cycleTime") << cycleTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("cycleCharge")
        << cycleCharge_ << token::END_STATEMENT << nl;
    os.writeKeyword("cycleEnergy")
        << cycleEnergy_ << token::END_STATEMENT << nl;
    os.writeKeyword("updatePeriods")
        << updatePeriods_ << token::END_STATEMENT << nl;
    os.writeKeyword("updateTime") << updateTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("updateEnergy")
        << updateEnergy_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rfPowerElectrode
    );
}

// ************************************************************************* //
