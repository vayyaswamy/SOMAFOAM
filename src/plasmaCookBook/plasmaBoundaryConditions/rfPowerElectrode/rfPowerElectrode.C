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
#include "Switch.H"

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
    powerTable_(0),
    powerRepeat_(false),
    tracking_("feedback"),
    waveform_("sine"),
    harmonics_(0),
    waveformTable_(0),
    amplitudeMax_(GREAT),
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
    updateTarget_(0),
    currentTarget_(-1),
    learned_(0),
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
    power_(dict.lookupOrDefault<scalar>("power", 0)),
    powerTable_(0),
    powerRepeat_(dict.lookupOrDefault<Switch>("powerRepeat", false)),
    tracking_(dict.lookupOrDefault<word>("tracking", "feedback")),
    waveform_(dict.lookupOrDefault<word>("waveform", "sine")),
    harmonics_(0),
    waveformTable_(0),
    amplitudeMax_(dict.lookupOrDefault<scalar>("amplitudeMax", GREAT)),
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
    updateTarget_(dict.lookupOrDefault<scalar>("updateTarget", 0)),
    currentTarget_(dict.lookupOrDefault<scalar>("currentTarget", -1)),
    learned_(0),
    curTimeIndex_(-1)
{
    if (dict.found("powerTable"))
    {
        dict.lookup("powerTable") >> powerTable_;
    }

    if (dict.found("harmonics"))
    {
        dict.lookup("harmonics") >> harmonics_;
    }

    if (dict.found("waveformTable"))
    {
        dict.lookup("waveformTable") >> waveformTable_;
    }

    if (dict.found("learned"))
    {
        dict.lookup("learned") >> learned_;
    }

    if
    (
        (waveform_ != "sine" && waveform_ != "harmonics" && waveform_ != "table")
     || (waveform_ == "harmonics" && harmonics_.empty())
     || (waveform_ == "table" && waveformTable_.size() < 2)
    )
    {
        FatalIOErrorIn
        (
            "rfPowerElectrode::rfPowerElectrode(...)",
            dict
        )   << "waveform must be sine, harmonics (with the entry harmonics)"
            << " or table (with the entry waveformTable)"
            << exit(FatalIOError);
    }

    if
    (
        (tracking_ != "feedback" && tracking_ != "learning")
     || (
            tracking_ == "learning"
         && (
                !powerRepeat_
             || powerTable_.size() < 2
             || mag(powerTable_[0].first()) > SMALL
            )
        )
     || (powerTable_.empty() && power_ <= 0)
    )
    {
        FatalIOErrorIn
        (
            "rfPowerElectrode::rfPowerElectrode(...)",
            dict
        )   << "Need power > 0 or a powerTable; tracking must be feedback"
            << " or learning, and learning needs a powerTable that starts"
            << " at time 0, with powerRepeat yes"
            << exit(FatalIOError);
    }

    if (tracking_ == "learning" && learned_.empty())
    {
        // One amplitude per period of the table, started from the
        // amplitude given for the largest power of the table
        const scalar tablePeriod = powerTable_[powerTable_.size() - 1].first();

        const label nPeriods = max(1, label(tablePeriod*frequency_ + 0.5));

        scalar maxPower = SMALL;

        forAll(powerTable_, i)
        {
            maxPower = max(maxPower, powerTable_[i].second());
        }

        learned_.setSize(nPeriods);

        forAll(learned_, i)
        {
            learned_[i] =
                amplitude_
               *Foam::sqrt
                (
                    max(targetPower((i + 0.5)/frequency_), scalar(0))/maxPower
                );
        }
    }

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
        frequency_ <= 0 || updateCycles_ < 1
     || relaxation_ <= 0 || relaxation_ > 1 || maxChange_ <= 0
     || capacitance_ <= 0
    )
    {
        FatalIOErrorIn
        (
            "rfPowerElectrode::rfPowerElectrode(...)",
            dict
        )   << "Need frequency > 0, updateCycles >= 1,"
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
    powerTable_(ptf.powerTable_),
    powerRepeat_(ptf.powerRepeat_),
    tracking_(ptf.tracking_),
    waveform_(ptf.waveform_),
    harmonics_(ptf.harmonics_),
    waveformTable_(ptf.waveformTable_),
    amplitudeMax_(ptf.amplitudeMax_),
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
    updateTarget_(ptf.updateTarget_),
    currentTarget_(ptf.currentTarget_),
    learned_(ptf.learned_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


Foam::rfPowerElectrode::rfPowerElectrode(const rfPowerElectrode& ptf)
:
    fixedValueFvPatchScalarField(ptf),
    frequency_(ptf.frequency_),
    power_(ptf.power_),
    powerTable_(ptf.powerTable_),
    powerRepeat_(ptf.powerRepeat_),
    tracking_(ptf.tracking_),
    waveform_(ptf.waveform_),
    harmonics_(ptf.harmonics_),
    waveformTable_(ptf.waveformTable_),
    amplitudeMax_(ptf.amplitudeMax_),
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
    updateTarget_(ptf.updateTarget_),
    currentTarget_(ptf.currentTarget_),
    learned_(ptf.learned_),
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
    powerTable_(ptf.powerTable_),
    powerRepeat_(ptf.powerRepeat_),
    tracking_(ptf.tracking_),
    waveform_(ptf.waveform_),
    harmonics_(ptf.harmonics_),
    waveformTable_(ptf.waveformTable_),
    amplitudeMax_(ptf.amplitudeMax_),
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
    updateTarget_(ptf.updateTarget_),
    currentTarget_(ptf.currentTarget_),
    learned_(ptf.learned_),
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
        currentTarget_ = targetPower((cycle_ + 0.5)/frequency_);

        if (tracking_ == "learning")
        {
            amplitude_ = learned_[cycle_ % learned_.size()];
        }

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

    // Target of the period that has ended and of the one that starts
    const scalar targetDone = targetPower((cycle_ + 0.5)/frequency_);
    const scalar targetNext = targetPower((cycleNow + 0.5)/frequency_);

    const label cycleDone = cycle_;
    cycle_ = cycleNow;

    const scalar meanCurrent = cycleCharge_/cycleTime_;
    const scalar meanPower = cycleEnergy_/cycleTime_;

    if (blockingCapacitor_ == "ideal")
    {
        // The mean current of the period charges the capacitor
        bias_ -= meanCurrent/frequency_/capacitance_;
    }

    if (tracking_ == "learning")
    {
        const label n = learned_.size();

        if (targetDone > 0)
        {
            scalar& a = learned_[cycleDone % n];
            a = min(a*correction(targetDone, meanPower), amplitudeMax_);
        }

        amplitude_ = learned_[cycleNow % n];
    }
    else
    {
        if (targetDone > 0)
        {
            updatePeriods_++;
            updateTime_ += cycleTime_;
            updateEnergy_ += cycleEnergy_;
            updateTarget_ += targetDone*cycleTime_;
        }

        if (updatePeriods_ >= updateCycles_)
        {
            amplitude_ *=
                correction
                (
                    updateTarget_/updateTime_,
                    updateEnergy_/updateTime_
                );

            updatePeriods_ = 0;
            updateTime_ = 0;
            updateEnergy_ = 0;
            updateTarget_ = 0;
        }

        // Feed-forward for a change of the target
        if (targetDone > 0 && targetNext > 0)
        {
            amplitude_ *= Foam::sqrt(targetNext/targetDone);
        }

        amplitude_ = min(amplitude_, amplitudeMax_);
    }

    currentTarget_ = targetNext;

    Info<< "rfPowerElectrode " << patch().name()
        << ": t " << tLast
        << " A " << (targetNext > 0 ? amplitude_ : 0)
        << " Vdc "
        << (blockingCapacitor_ == "physical" ? bias_ - charge_/capacitance_ : bias_)
        << " P " << meanPower
        << " I " << meanCurrent
        << " target " << targetDone << endl;

    cycleTime_ = 0;
    cycleCharge_ = 0;
    cycleEnergy_ = 0;
}


Foam::scalar Foam::rfPowerElectrode::correction
(
    const scalar target,
    const scalar measured
) const
{
    scalar factor = 1 + maxChange_;

    if (measured > 0)
    {
        factor = pow(target/measured, 0.5*relaxation_);
        factor = max(min(factor, 1 + maxChange_), 1/(1 + maxChange_));
    }

    return factor;
}


Foam::scalar Foam::rfPowerElectrode::interpolate
(
    const List<Tuple2<scalar, scalar> >& table,
    const scalar x
)
{
    const label n = table.size();

    if (x <= table[0].first())
    {
        return table[0].second();
    }

    for (label i = 1; i < n; i++)
    {
        if (x <= table[i].first())
        {
            const scalar dx = table[i].first() - table[i - 1].first();

            if (dx <= VSMALL)
            {
                return table[i].second();
            }

            const scalar w = (x - table[i - 1].first())/dx;

            return (1 - w)*table[i - 1].second() + w*table[i].second();
        }
    }

    return table[n - 1].second();
}


Foam::scalar Foam::rfPowerElectrode::targetPower(const scalar t) const
{
    if (powerTable_.empty())
    {
        return power_;
    }

    scalar tt = t;

    if (powerRepeat_)
    {
        const scalar tablePeriod = powerTable_[powerTable_.size() - 1].first();

        tt = t - tablePeriod*floor(t/tablePeriod);
    }

    return interpolate(powerTable_, tt);
}


Foam::scalar Foam::rfPowerElectrode::shape(const scalar t) const
{
    const scalar twoPi = 2*mathematicalConstant::pi;

    if (waveform_ == "sine")
    {
        return Foam::sin(twoPi*frequency_*t);
    }
    else if (waveform_ == "harmonics")
    {
        scalar g = 0;

        forAll(harmonics_, i)
        {
            const vector& h = harmonics_[i];

            g += h.y()*Foam::sin(twoPi*h.x()*frequency_*t + h.z()*twoPi/360.0);
        }

        return g;
    }

    const scalar x = frequency_*t;

    return interpolate(waveformTable_, x - floor(x));
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

    if (currentTarget_ < 0)
    {
        // Restart from a field written without it
        currentTarget_ = targetPower(runTime.value());
    }

    // No voltage in periods with zero target power
    lastVoltage_ =
        (currentTarget_ > 0 ? amplitude_*shape(runTime.value()) : 0) + Vdc;

    operator==(lastVoltage_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::rfPowerElectrode::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    if (powerTable_.empty())
    {
        os.writeKeyword("power") << power_ << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("powerTable")
            << powerTable_ << token::END_STATEMENT << nl;
        os.writeKeyword("powerRepeat")
            << Switch(powerRepeat_) << token::END_STATEMENT << nl;
    }

    os.writeKeyword("tracking") << tracking_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveform") << waveform_ << token::END_STATEMENT << nl;

    if (waveform_ == "harmonics")
    {
        os.writeKeyword("harmonics")
            << harmonics_ << token::END_STATEMENT << nl;
    }
    else if (waveform_ == "table")
    {
        os.writeKeyword("waveformTable")
            << waveformTable_ << token::END_STATEMENT << nl;
    }

    if (amplitudeMax_ < GREAT)
    {
        os.writeKeyword("amplitudeMax")
            << amplitudeMax_ << token::END_STATEMENT << nl;
    }
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
    os.writeKeyword("updateTarget")
        << updateTarget_ << token::END_STATEMENT << nl;
    os.writeKeyword("currentTarget")
        << currentTarget_ << token::END_STATEMENT << nl;

    if (learned_.size())
    {
        os.writeKeyword("learned") << learned_ << token::END_STATEMENT << nl;
    }

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
