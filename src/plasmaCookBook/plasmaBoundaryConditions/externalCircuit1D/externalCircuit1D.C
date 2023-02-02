/*---------------------------------------------------------------------------*\
Copyright (C) 2023 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "externalCircuit1D.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mathematicalConstants.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from patch and internal field
Foam::externalCircuit1D::externalCircuit1D
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
	distance_(0.0),
	init_(false),
	oldValue_(0.0),
	newValue_(0.0)
	{

	}

// Construct from patch, internal field and dictionary
Foam::externalCircuit1D::externalCircuit1D
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
	distance_(0.0), 
	init_(false),
	oldValue_(0.0),
	newValue_(0.0)
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

		// calculate mesh length, for a 1D mesh in the x direction 
		const faceList & ff = this->patch().boundaryMesh().mesh().faces(); 
		const pointField & pp = this->patch().boundaryMesh().mesh().points(); 

		forAll(this->patch().boundaryMesh().mesh().C(), celli)
		{
			const cell & cc = this->patch().boundaryMesh().mesh().cells()[celli]; 
			labelList pLabels(cc.labels(ff)); 
			pointField pLocal(pLabels.size(), vector::zero); 
			forAll (pLabels, pointi)
			    pLocal[pointi] = pp[pLabels[pointi]]; 

			scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0)); // And similar for yDim and zDim 
			distance_ += xDim;
		}	
	}

// Construct by mapping given externalCircuit1D onto a new patch
Foam::externalCircuit1D::externalCircuit1D
(
    const externalCircuit1D& ptf,
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
    distance_(ptf.distance_),
	init_(ptf.init_),
	oldValue_(ptf.oldValue_),
	newValue_(ptf.newValue_)
    {

    }


// Construct as copy
Foam::externalCircuit1D::externalCircuit1D
(
	const externalCircuit1D& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_),
    R_(tppsf.R_),
    distance_(tppsf.distance_),
	init_(tppsf.init_),
	oldValue_(tppsf.oldValue_),
	newValue_(tppsf.newValue_)
    {

    }

// Construct a copy setting internal field reference
Foam::externalCircuit1D::externalCircuit1D
(
    const externalCircuit1D& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    amplitude_(tppsf.amplitude_),
    modelName_(tppsf.modelName_),
    frequency_(tppsf.frequency_),
    bias_(tppsf.bias_),
    R_(tppsf.R_),
    distance_(tppsf.distance_),
    init_(tppsf.init_),
    oldValue_(tppsf.oldValue_),
    newValue_(tppsf.newValue_)
    {

    }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCircuit1D::updateCoeffs()
{
	if (updated())
	{
		return;
	}

	if(init_ == true)
	{
		// Info << "boundary coundition true: " << oldValue_ << endl;
		operator==(oldValue_);
	}
	else if(init_ == false)
	{
		oldValue_=0.75*amplitude_[0];
		operator==(oldValue_);
		// Info << "initializing boundary condition: " << oldValue_ << endl;
		init_ = true;
	}

	// obtain area of operation
	// patch() ->  fvPatchField.H
	// boundaryMesh() -> fvPatch.H
	// mesh() -> fvBoundaryMesh.H
	// boundaryMesh() -> polyMesh.H
	// findPatchID() -> polyBoundaryMesh.H  
	label electrode = this->patch().boundaryMesh().mesh().boundaryMesh().findPatchID("electrode");
	scalar areaMesh = this->patch().boundaryMesh().mesh().magSf().boundaryField()[electrode][0];

	scalar resistance = areaMesh*R_;
	scalar epsiloni = 0;

	// Following code block is to generate the future boundary value 
	if(modelName_ == "cosFrequencyModulated")
	{
		scalar cosineFunction = amplitude_[0]*Foam::cos(2*mathematicalConstant::pi*frequency_*this->db().time().value());

		const volVectorField& currentDensity = db().objectRegistry::lookupObject<volVectorField>("conductionCurrent");
		const volScalarField& epsilon = db().objectRegistry::lookupObject<volScalarField>("epsiloni");
		const scalarField& x_currentDensity = currentDensity.component(0);

		scalar netCurrent = 0; 
		scalar sumCurrent = 0;
		const faceList & ff = this->patch().boundaryMesh().mesh().faces(); 
		const pointField & pp = this->patch().boundaryMesh().mesh().points();

		forAll(this->patch().boundaryMesh().mesh().C(), celli)
		{
			const cell & cc = this->patch().boundaryMesh().mesh().cells()[celli]; 
			labelList pLabels(cc.labels(ff)); 
			pointField pLocal(pLabels.size(), vector::zero); 
			forAll (pLabels, pointi)
			    pLocal[pointi] = pp[pLabels[pointi]]; 

			scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0)); // And similar for yDim and zDim 
			sumCurrent += x_currentDensity[celli]*xDim;
			epsiloni += epsilon[celli]*xDim;
		}

		epsiloni = (1/distance_)*epsiloni;

		scalar coeffs = (distance_*this->db().time().deltaTValue())/(epsiloni*resistance);

		scalar dropResistance = (cosineFunction-oldValue_)/(1+coeffs);

		netCurrent = (1/distance_)*sumCurrent;

		newValue_ = oldValue_+(dropResistance-resistance*netCurrent)*coeffs;
		oldValue_ = newValue_;
	}
	else
	{
        FatalErrorIn
        (
            "externalCircuit1D::updateCoeffs()"
        )   << " model name inconsitent, model = " << modelName_
            << exit(FatalError);
	}

	fixedValueFvPatchScalarField::updateCoeffs();
}

// write boundary conditions parameters when writeControl requires so 
void Foam::externalCircuit1D::write(Ostream& os) const
{
	fvPatchScalarField::write(os);

	os.writeKeyword("model") << modelName_ << token::END_STATEMENT << nl;
	amplitude_.writeEntry("amplitude", os);
	os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
	os.writeKeyword("bias") << bias_ << token::END_STATEMENT << nl;
	os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;

	writeEntry("value", os);
}

// macro required for compilation purporses 
// fvPatchField.H
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, externalCircuit1D);
}
