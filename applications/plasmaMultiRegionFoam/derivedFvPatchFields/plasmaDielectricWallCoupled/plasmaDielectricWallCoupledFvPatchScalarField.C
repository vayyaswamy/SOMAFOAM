/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaDielectricWallCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaDielectricWallCoupledFvPatchScalarField::
plasmaDielectricWallCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    neighbourPatchName_("undefined-neighbourPatchName")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


Foam::plasmaDielectricWallCoupledFvPatchScalarField::
plasmaDielectricWallCoupledFvPatchScalarField
(
    const plasmaDielectricWallCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourPatchName_(ptf.neighbourPatchName_)
{}


Foam::plasmaDielectricWallCoupledFvPatchScalarField::
plasmaDielectricWallCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    neighbourPatchName_(dict.lookup("neighbourPatchName"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "plasmaDielectricWallCoupledFvPatchScalarField::"
            "plasmaDielectricWallCoupledFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


Foam::plasmaDielectricWallCoupledFvPatchScalarField::
plasmaDielectricWallCoupledFvPatchScalarField
(
    const plasmaDielectricWallCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    neighbourPatchName_(wtcsf.neighbourPatchName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvPatchScalarField&
Foam::plasmaDielectricWallCoupledFvPatchScalarField::sigma() const
{
    const fvPatchField<scalar>& tsurfC =
        patch().lookupPatchField<volScalarField, scalar>("surfC");

    return tsurfC;
}

//const Foam::fvPatchVectorField&
//Foam::plasmaDielectricWallCoupledFvPatchScalarField::E() const
//{
//    const fvPatchField<vector>& tE =
//        patch().lookupPatchField<volVectorField, vector>("E");

//    return tE;
//}

const Foam::fvPatchScalarField&
Foam::plasmaDielectricWallCoupledFvPatchScalarField::epsilon() const
{
    const fvPatchField<scalar>& tEps =
        patch().lookupPatchField<volScalarField, scalar>("epsilon");

    return tEps;
}

void Foam::plasmaDielectricWallCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the directMappedPatchBase
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();

    tmp<scalarField> intFld = patchInternalField();

    const plasmaDielectricWallCoupledFvPatchScalarField& nbrField =
    refCast<const plasmaDielectricWallCoupledFvPatchScalarField>
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            neighbourFieldName_
        )
    );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld = nbrField.patchInternalField();
    mapDistribute::distribute
    (
        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrIntFld
    );

    scalarField nbrSigma = nbrField.sigma();
    mapDistribute::distribute
    (
        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrSigma
    );

    scalarField nbrEpsilon = nbrField.epsilon()*nbrPatch.deltaCoeffs();
    mapDistribute::distribute
    (
        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrEpsilon
    );

    tmp<scalarField> myEpsilon = epsilon()*patch().deltaCoeffs();

	if (neighbourPatchName_ == "plasma")
	{
		this->refValue() = nbrIntFld;

		this->refGrad() = - sigma()/epsilon();

		this->valueFraction() = nbrEpsilon / (nbrEpsilon + myEpsilon());
	}
	else if (neighbourPatchName_ == "dielectric")
	{
		this->refValue() = nbrIntFld;

		this->refGrad() = - nbrSigma/epsilon();

		this->valueFraction() = nbrEpsilon / (nbrEpsilon + myEpsilon());
	}


//	if (neighbourPatchName_ == "plasma")
//	{
//		this->refValue() = 0.0;

//		this->refGrad() = (sigma()-nbrEpsilon*(nbrElectricField & patch().nf()))/epsilon();

//		this->valueFraction() = 0.0;
//	}
//	else if (neighbourPatchName_ == "dielectric")
//	{
//		this->refValue() = nbrIntFld;

//		this->refGrad() = 0.0;

//		this->valueFraction() = 1.0;
//	}

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::plasmaDielectricWallCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName")<< neighbourPatchName_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    plasmaDielectricWallCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
