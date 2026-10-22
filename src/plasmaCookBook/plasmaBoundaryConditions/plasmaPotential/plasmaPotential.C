/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

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
#include "DynamicList.H"

#include <fstream>
#include <sstream>

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
    bias_(0.0),
    frequencyTable_(0),
    envelopeTable_(0),
    phase_(0.0),
    voltageTable_(0),
    voltageTableFile_(),
    dtControl_(false),
    maxRelChange_(0.05),
    currentFloor_(0),
    maxIncrease_(1.2),
    sampleInterval_(0),
    lastTimeIndex_(-1),
    nSamples_(0),
    tA_(0),
    tB_(0),
    qB_(0),
    iB_(0),
    current_(0),
    targetDeltaT_(GREAT),
    suggestedDeltaT_(GREAT)
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
    bias_(dict.lookupOrDefault<scalar>("bias", 0.0)),
    frequencyTable_(0),
    envelopeTable_(0),
    phase_(dict.lookupOrDefault<scalar>("phase", 0.0)),
    voltageTable_(0),
    voltageTableFile_(),
    dtControl_(false),
    maxRelChange_(0.05),
    currentFloor_(0),
    maxIncrease_(1.2),
    sampleInterval_(0),
    lastTimeIndex_(-1),
    nSamples_(0),
    tA_(0),
    tB_(0),
    qB_(0),
    iB_(0),
    current_(0),
    targetDeltaT_(GREAT),
    suggestedDeltaT_(GREAT)
{
    readModel(dict);
    readTimeStepControl(dict);

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
    bias_(ptf.bias_),
    frequencyTable_(ptf.frequencyTable_),
    envelopeTable_(ptf.envelopeTable_),
    phase_(ptf.phase_),
    voltageTable_(ptf.voltageTable_),
    voltageTableFile_(ptf.voltageTableFile_),
    dtControl_(ptf.dtControl_),
    maxRelChange_(ptf.maxRelChange_),
    currentFloor_(ptf.currentFloor_),
    maxIncrease_(ptf.maxIncrease_),
    sampleInterval_(ptf.sampleInterval_),
    lastTimeIndex_(ptf.lastTimeIndex_),
    nSamples_(ptf.nSamples_),
    tA_(ptf.tA_),
    tB_(ptf.tB_),
    qB_(ptf.qB_),
    iB_(ptf.iB_),
    current_(ptf.current_),
    targetDeltaT_(ptf.targetDeltaT_),
    suggestedDeltaT_(ptf.suggestedDeltaT_)
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
    bias_(tppsf.bias_),
    frequencyTable_(tppsf.frequencyTable_),
    envelopeTable_(tppsf.envelopeTable_),
    phase_(tppsf.phase_),
    voltageTable_(tppsf.voltageTable_),
    voltageTableFile_(tppsf.voltageTableFile_),
    dtControl_(tppsf.dtControl_),
    maxRelChange_(tppsf.maxRelChange_),
    currentFloor_(tppsf.currentFloor_),
    maxIncrease_(tppsf.maxIncrease_),
    sampleInterval_(tppsf.sampleInterval_),
    lastTimeIndex_(tppsf.lastTimeIndex_),
    nSamples_(tppsf.nSamples_),
    tA_(tppsf.tA_),
    tB_(tppsf.tB_),
    qB_(tppsf.qB_),
    iB_(tppsf.iB_),
    current_(tppsf.current_),
    targetDeltaT_(tppsf.targetDeltaT_),
    suggestedDeltaT_(tppsf.suggestedDeltaT_)
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
    bias_(tppsf.bias_),
    frequencyTable_(tppsf.frequencyTable_),
    envelopeTable_(tppsf.envelopeTable_),
    phase_(tppsf.phase_),
    voltageTable_(tppsf.voltageTable_),
    voltageTableFile_(tppsf.voltageTableFile_),
    dtControl_(tppsf.dtControl_),
    maxRelChange_(tppsf.maxRelChange_),
    currentFloor_(tppsf.currentFloor_),
    maxIncrease_(tppsf.maxIncrease_),
    sampleInterval_(tppsf.sampleInterval_),
    lastTimeIndex_(tppsf.lastTimeIndex_),
    nSamples_(tppsf.nSamples_),
    tA_(tppsf.tA_),
    tB_(tppsf.tB_),
    qB_(tppsf.qB_),
    iB_(tppsf.iB_),
    current_(tppsf.current_),
    targetDeltaT_(tppsf.targetDeltaT_),
    suggestedDeltaT_(tppsf.suggestedDeltaT_)
{}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plasmaPotential::readModel(const dictionary& dict)
{
    if (modelName_ == "chirp")
    {
        dict.lookup("frequencyTable") >> frequencyTable_;

        if (dict.found("envelopeTable"))
        {
            dict.lookup("envelopeTable") >> envelopeTable_;
        }

        if (frequencyTable_.empty() || frequencyTable_[0].first() < 0)
        {
            FatalIOErrorIn("plasmaPotential::readModel", dict)
                << "chirp: frequencyTable must be non-empty and start at"
                << " t >= 0" << exit(FatalIOError);
        }
    }
    else if (modelName_ == "table")
    {
        if (dict.found("voltageTableFile"))
        {
            voltageTableFile_ = fileName(dict.lookup("voltageTableFile"));

            fileName f(voltageTableFile_);
            f.expand();

            if (f.empty() || f[0] != '/')
            {
                // Relative to the case directory (also in parallel)
                const Time& runTime = this->db().time();
                f = runTime.processorCase()
                  ? runTime.path()/".."/f
                  : runTime.path()/f;
            }

            std::ifstream is(f.c_str());

            if (!is.good())
            {
                FatalIOErrorIn("plasmaPotential::readModel", dict)
                    << "table: cannot open voltageTableFile " << f
                    << exit(FatalIOError);
            }

            DynamicList<Tuple2<scalar, scalar> > rows;
            std::string line;

            while (std::getline(is, line))
            {
                std::istringstream ls(line);
                double t, v;

                // Skip blank lines and lines that do not start with two
                // numbers (comments, headers)
                if (ls >> t >> v)
                {
                    rows.append(Tuple2<scalar, scalar>(t, v));
                }
            }

            voltageTable_.setSize(rows.size());

            forAll(rows, i)
            {
                voltageTable_[i] = rows[i];
            }
        }
        else
        {
            dict.lookup("voltageTable") >> voltageTable_;
        }

        if (voltageTable_.size() < 2)
        {
            FatalIOErrorIn("plasmaPotential::readModel", dict)
                << "table: give voltageTable or voltageTableFile with at"
                << " least two entries" << exit(FatalIOError);
        }
    }
    else if
    (
        modelName_ != "directCurrent"
     && modelName_ != "sinFrequencyModulated"
     && modelName_ != "cosFrequencyModulated"
    )
    {
        FatalIOErrorIn("plasmaPotential::readModel", dict)
            << "unknown model " << modelName_ << ". Valid models:"
            << " directCurrent, sinFrequencyModulated, cosFrequencyModulated,"
            << " chirp, table" << exit(FatalIOError);
    }

    // Times must not decrease (equal times give a step)
    const List<Tuple2<scalar, scalar> >* tables[3] =
        {&frequencyTable_, &envelopeTable_, &voltageTable_};

    for (label k = 0; k < 3; k++)
    {
        const List<Tuple2<scalar, scalar> >& tbl = *tables[k];

        for (label i = 1; i < tbl.size(); i++)
        {
            if (tbl[i].first() < tbl[i - 1].first())
            {
                FatalIOErrorIn("plasmaPotential::readModel", dict)
                    << "table times must not decrease (entry " << i << ")"
                    << exit(FatalIOError);
            }
        }
    }
}


Foam::scalar Foam::plasmaPotential::interpolate
(
    const List<Tuple2<scalar, scalar> >& table,
    const scalar t
)
{
    const label n = table.size();

    if (t <= table[0].first())
    {
        return table[0].second();
    }

    if (t >= table[n - 1].first())
    {
        return table[n - 1].second();
    }

    // Bisection: table[lo].first() <= t < table[hi].first()
    label lo = 0;
    label hi = n - 1;

    while (hi - lo > 1)
    {
        const label mid = (lo + hi)/2;

        if (table[mid].first() <= t)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    const scalar dt = table[hi].first() - table[lo].first();

    if (dt <= VSMALL)
    {
        return table[hi].second();
    }

    const scalar w = (t - table[lo].first())/dt;

    return (1 - w)*table[lo].second() + w*table[hi].second();
}


Foam::scalar Foam::plasmaPotential::chirpPhase(const scalar t) const
{
    // Number of cycles from t = 0: the integral of the piecewise-linear f(t),
    // with f held at its end values outside the table
    const List<Tuple2<scalar, scalar> >& tbl = frequencyTable_;
    const label n = tbl.size();

    scalar cycles = 0;

    if (t <= tbl[0].first())
    {
        cycles = tbl[0].second()*t;
    }
    else
    {
        cycles = tbl[0].second()*tbl[0].first();

        label i = 1;

        for (; i < n; i++)
        {
            const scalar ta = tbl[i - 1].first();
            const scalar tb = tbl[i].first();
            const scalar fa = tbl[i - 1].second();
            const scalar fb = tbl[i].second();

            if (tb - ta <= VSMALL)
            {
                continue;
            }

            if (t <= tb)
            {
                const scalar s = t - ta;
                cycles += fa*s + 0.5*(fb - fa)/(tb - ta)*s*s;
                break;
            }

            cycles += 0.5*(fa + fb)*(tb - ta);
        }

        if (i == n)
        {
            cycles += tbl[n - 1].second()*(t - tbl[n - 1].first());
        }
    }

    return
        2*mathematicalConstant::pi*cycles
      + phase_*mathematicalConstant::pi/180.0;
}


void Foam::plasmaPotential::readTimeStepControl(const dictionary& dict)
{
    if (!dict.found("timeStepControl"))
    {
        return;
    }

    const dictionary& d = dict.subDict("timeStepControl");

    dtControl_ = d.lookupOrDefault<Switch>("active", true);
    maxRelChange_ = d.lookupOrDefault<scalar>("maxRelChange", 0.05);
    currentFloor_ = d.lookupOrDefault<scalar>("currentFloor", 0.0);
    maxIncrease_ = d.lookupOrDefault<scalar>("maxIncrease", 1.2);
    sampleInterval_ = d.lookupOrDefault<scalar>("sampleInterval", 0.0);

    if
    (
        maxRelChange_ <= 0 || currentFloor_ < 0 || maxIncrease_ < 1
     || sampleInterval_ < 0
    )
    {
        FatalIOErrorIn("plasmaPotential::readTimeStepControl", d)
            << "timeStepControl needs maxRelChange > 0, currentFloor >= 0,"
            << " maxIncrease >= 1 and sampleInterval >= 0"
            << exit(FatalIOError);
    }
}


Foam::scalar Foam::plasmaPotential::patchCharge() const
{
    // eps on the patch: the region's epsilon field, or vacuum if absent
    scalarField epsp(size(), 8.85418782e-12);

    if (db().foundObject<volScalarField>("epsilon"))
    {
        epsp = patch().lookupPatchField<volScalarField, scalar>("epsilon");
    }

    const scalarField dPhidn
    (
        patch().deltaCoeffs()*(scalarField(*this) - patchInternalField())
    );

    return gSum(epsp*dPhidn*patch().magSf());
}


void Foam::plasmaPotential::sampleCurrent()
{
    const Time& runTime = db().time();

    if (!dtControl_ || runTime.timeIndex() == lastTimeIndex_)
    {
        return;
    }

    lastTimeIndex_ = runTime.timeIndex();

    // The first update in a time step comes before this step's solution, so
    // the patch and cell values still belong to the end of the previous step
    const scalar dt = runTime.deltaTValue();
    const scalar t = runTime.value() - dt;

    if (nSamples_ == 0)
    {
        // First sample (start or restart): hold the step until the current
        // has been measured over two windows
        tA_ = t;
        tB_ = t;
        qB_ = patchCharge();
        targetDeltaT_ = dt;
        nSamples_ = 1;
    }
    else if (t - tB_ >= max(sampleInterval_, SMALL*mag(t)))
    {
        // End of a sampling window: current over it, and its rate of change
        // between the centres of this window and the previous one
        const scalar q = patchCharge();
        const scalar I = (q - qB_)/(t - tB_);

        if (nSamples_ >= 2)
        {
            const scalar rate = mag(I - iB_)/(0.5*(t - tA_));

            targetDeltaT_ =
                maxRelChange_*(mag(I) + currentFloor_)/(rate + VSMALL);
        }

        tA_ = tB_;
        tB_ = t;
        qB_ = q;
        iB_ = I;
        current_ = I;
        nSamples_++;
    }

    // Drop to the target at once; grow by at most maxIncrease per step
    suggestedDeltaT_ = min(targetDeltaT_, maxIncrease_*dt);
}


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

    sampleCurrent();

    if (modelName_ == "directCurrent")
    {
        operator==(amplitude_);
    }
    else if (modelName_ == "cosFrequencyModulated")
    {
        operator==(amplitude_*Foam::cos(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_);
    }
    else if (modelName_ == "sinFrequencyModulated")
    {
        operator==(amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*this->db().time().value()) + bias_);
    }
    else if (modelName_ == "chirp")
    {
        const scalar t = this->db().time().value();
        const scalar e =
            envelopeTable_.size() ? interpolate(envelopeTable_, t) : 1.0;

        operator==(amplitude_*e*Foam::sin(chirpPhase(t)) + bias_);
    }
    else if (modelName_ == "table")
    {
        const scalar t = this->db().time().value();

        operator==(amplitude_*interpolate(voltageTable_, t) + bias_);
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

    if (modelName_ == "chirp")
    {
        os.writeKeyword("frequencyTable")
            << frequencyTable_ << token::END_STATEMENT << nl;

        if (envelopeTable_.size())
        {
            os.writeKeyword("envelopeTable")
                << envelopeTable_ << token::END_STATEMENT << nl;
        }

        os.writeKeyword("phase")
            << phase_ << token::END_STATEMENT << nl;
    }
    else if (modelName_ == "table")
    {
        // Keep the file name rather than the data, so that restarts read
        // the same file and the time directories stay small
        if (voltageTableFile_.size())
        {
            os.writeKeyword("voltageTableFile")
                << voltageTableFile_ << token::END_STATEMENT << nl;
        }
        else
        {
            os.writeKeyword("voltageTable")
                << voltageTable_ << token::END_STATEMENT << nl;
        }
    }

    if (dtControl_)
    {
        os << indent << "timeStepControl" << nl
           << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("active")
            << dtControl_ << token::END_STATEMENT << nl;
        os.writeKeyword("maxRelChange")
            << maxRelChange_ << token::END_STATEMENT << nl;
        os.writeKeyword("currentFloor")
            << currentFloor_ << token::END_STATEMENT << nl;
        os.writeKeyword("maxIncrease")
            << maxIncrease_ << token::END_STATEMENT << nl;
        os.writeKeyword("sampleInterval")
            << sampleInterval_ << token::END_STATEMENT << nl;
        os << decrIndent << indent << token::END_BLOCK << nl;
    }

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, plasmaPotential);
}
// ************************************************************************* //
