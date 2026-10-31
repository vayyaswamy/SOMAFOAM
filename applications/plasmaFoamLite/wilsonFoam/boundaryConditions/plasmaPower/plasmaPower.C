/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors
(This code was written by Abhishek Kumar Verma a former PhD student at UC Merced)
License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaPower.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mathematicalConstants.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::plasmaPower::
plasmaPower
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    //amplitude_(p.size(), 0.0),
    amplitude_(0.0),
    //modelName_("directCurrent"),
    frequency_(0.0),
    power_(0.0),
    nCycles_(1.0),
    dampingFactor_(0.0),
    bias_(0.0)
{}


Foam::plasmaPower::
plasmaPower
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    //amplitude_("amplitude", dict, p.size()),
    amplitude_(dict.lookupOrDefault<scalar>("amplitude", 0.0)),
    //modelName_(dict.lookupOrDefault<word>("model", "directCurrent")),
    frequency_(dict.lookupOrDefault<scalar>("frequency", 0.0)),
    power_(dict.lookupOrDefault<scalar>("power", 0.0)),
    nCycles_(dict.lookupOrDefault<scalar>("nCycles", 0.0)),
    dampingFactor_(dict.lookupOrDefault<scalar>("dampingFactor", 0.0)),
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

Foam::plasmaPower::
plasmaPower
(
    const plasmaPower& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    //amplitude_(ptf.amplitude_, mapper),
    amplitude_(ptf.amplitude_),
    //modelName_(ptf.modelName_),
    frequency_(ptf.frequency_),
    power_(ptf.power_),
    nCycles_(ptf.nCycles_),
    dampingFactor_(ptf.dampingFactor_),
    bias_(ptf.bias_)
{}


Foam::plasmaPower::
plasmaPower
(
    const plasmaPower& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    //amplitude_(tppsf.amplitude_),
    amplitude_(tppsf.amplitude_),
    //modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    power_(tppsf.power_),
    nCycles_(tppsf.nCycles_),
    dampingFactor_(tppsf.dampingFactor_),
    bias_(tppsf.bias_)
{}


Foam::plasmaPower::
plasmaPower
(
    const plasmaPower& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    //amplitude_(tppsf.amplitude_),
    amplitude_(tppsf.amplitude_),
    //modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    power_(tppsf.power_),
    nCycles_(tppsf.nCycles_),
    dampingFactor_(tppsf.dampingFactor_),
    bias_(tppsf.bias_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaPower::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}

void Foam::plasmaPower::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}
//=================================================




void Foam::plasmaPower::updateCoeffs()
{
    //get power value from field (maybe not good method but works perfectly)
    const fvPatchField<scalar>& powerR=
        patch().lookupPatchField<volScalarField, scalar>("power");
    scalar currentPower = gAverage(powerR);
    sumPowerr_ = sumPowerr_ + currentPower;
    //--what is current time
    currentTimee = this->db().time().value(); 
    //Info << this->db().time().value() << endl;
    //--what was previous time
    Info << "previous time" << previousTimee << endl;
    Info << "current time" << currentTimee << endl;
    Info << "difference in time (curr-prev)" <<(currentTimee - previousTimee) << "cyccles.freq" << nCycles_/frequency_ << endl;
    Info << "------------" << nCycles_ << endl; 
    Info << "Damping Factor" <<dampingFactor_ << endl;
    //--checking if the cycle has crossed
    nCountter =nCountter+1;
    Info << "count ========" <<nCountter << endl;

    string filename("powerVoltageLogFile");
    Info << "========= sum power===" << sumPowerr_/nCountter << " power" << power_ << endl;
    Info << "====== amplitude =" << amplitude_ << endl;
    Info << "======= times correct" <<timesCorrect << endl;
    if ((currentTimee - previousTimee)> nCycles_/frequency_){
        timesCorrect = timesCorrect+1;
        Info << "la bhayo" << endl;
        previousTimee = currentTimee;
        averagePower = sumPowerr_/nCountter;
        Info << "la bhayo" << endl;
        nCountter= 0;
        sumPowerr_ = 0;
        //changing the amplitude
        scalar oldAmplitude = amplitude_;
        
        amplitude_ = oldAmplitude*(1.0-dampingFactor_*((averagePower/power_)-1.0));
        Info << "la bhayo3" << endl;
        if((amplitude_/oldAmplitude) >= 1.1)
        {
            amplitude_ = 1.1*oldAmplitude;
        }
        else if((amplitude_/oldAmplitude) <= 0.9)
        {
            amplitude_ = 0.9*oldAmplitude;
        }
        Info << "power ko antya samma aayo" << endl;
        //write -----------------------------------
        ofstream myfile;
        myfile.open ("power_voltage_log.txt",std::ios::app);
        myfile << "time: " << currentTimee << "\t\tAmplitude: " << amplitude_ <<  "\t\tPower: " << averagePower << "\n" << endl;
        myfile.close();
    			//if (Pstream::master())
		//	{
		/*	   OFstream& powerLogFile = *powerLogFilePtr_;
			   int width = 20;
			   powerLogFile << currentTimee;
			   powerLogFile.width(width);
			   powerLogFile << averagePower;
			   powerLogFile.width(width);
			   powerLogFile << amplitude_;
			   powerLogFile << endl;
		*///	}
        // --------------------------------------
    }
    if (updated())
    {
        return;
    }
    operator==(amplitude_*Foam::cos(2*22/7*frequency_*this->db().time().value()) );
    Info << "voltageValue" << (amplitude_*Foam::sin(2*22/7*frequency_*this->db().time().value()) ) << endl;
    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::plasmaPower::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    //os.writeKeyword("model")
    //    << modelName_ << token::END_STATEMENT << nl;
    //amplitude_.writeEntry("amplitude", os);
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("nCycles")
        << nCycles_ << token::END_STATEMENT << nl;
    os.writeKeyword("power")
        << power_ << token::END_STATEMENT << nl;
    os.writeKeyword("dampingFactor")
        << dampingFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("bias")
        << bias_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, plasmaPower);
}
// ************************************************************************* //
