/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors
License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.
\*---------------------------------------------------------------------------*/
#include "mappingPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
//****Chenhui edit start
//****Note: directMappedPatchBase has been dropped by OpenFoam
//#include "directMappedPatchBase.H"
#include "mappedPatchBase.H"
//****Chenhui edit end

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::mappingPotential::
mappingPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    neighbourPatchName_("undefined-neighbourPatchName"),
    neighbourPatchType_("undefined-neighbourPatchType")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

Foam::mappingPotential::
mappingPotential
(
    const mappingPotential& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourPatchName_(ptf.neighbourPatchName_),
    neighbourPatchType_(ptf.neighbourPatchType_)
{}

Foam::mappingPotential::
mappingPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    neighbourPatchName_(dict.lookup("neighbourPatchName")),
    neighbourPatchType_(dict.lookup("neighbourPatchType"))
{
//****Chenhui edit start
//****Note: directMappedPatchBase has been dropped from OpenFOAM
//    if (!isA<directMappedPatchBase>(this->patch().patch()))
    if (!isA<mappedPatchBase>(this->patch().patch()))
//***Chenhui edit end
    {
        FatalErrorIn
        (
            "mappingPotential::"
            "mappingPotential\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
//****Chenhui edit start
//****Note: directMappedPatchBase has been dropped from OpenFOAM
//            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "' not type '" << mappedPatchBase::typeName << "'"
//****Chenhui edit end
            << "\n    for patch " << p.name()
//***Chenhui edit start
//***Note: dimensionedInternalField is not defined in this scope
//            << " of field " << dimensionedInternalField().name()
//            << " in file " << dimensionedInternalField().objectPath()
//***Chenhui edit end
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    if (dict.found("refValue"))
    {
        Info << "refBhetayo" << endl;
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        Info << "refBhetayo===========" << endl;
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}

Foam::mappingPotential::
mappingPotential
(
    const mappingPotential& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    neighbourPatchName_(wtcsf.neighbourPatchName_),
    neighbourPatchType_(wtcsf.neighbourPatchType_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvPatchScalarField&
Foam::mappingPotential::sigma() const
{
    const fvPatchField<scalar>& tsurfC =
        patch().lookupPatchField<volScalarField, scalar>("surfC");

    return tsurfC;
}

const Foam::fvPatchVectorField&
Foam::mappingPotential::EField() const
{
    const fvPatchField<vector>& tE =
        patch().lookupPatchField<volVectorField, vector>("EField");
    return tE;
}

const Foam::fvPatchScalarField&
Foam::mappingPotential::epsilon() const
{
    const fvPatchField<scalar>& tEps =
        patch().lookupPatchField<volScalarField, scalar>("epsilon");
    return tEps;
}

const Foam::fvPatchScalarField&
Foam::mappingPotential::EPotential() const
{
    const fvPatchField<scalar>& tEPotential =
        patch().lookupPatchField<volScalarField, scalar>("EPotential");
    return tEPotential;
}


void Foam::mappingPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }
//****Chenhui edit start
//****Note: directMappedPatchBase has been dropped from OpenFOAM
//    // Get the coupling information from the directMappedPatchBase
//    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
//***Chenhui edit end
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

//***Chenhui edit start
//***Note: distMap is not used in OpenFOAM version of the BC
//    // Force recalculation of mapping and schedule
//    const mapDistribute& distMap = mpp.map();
//***Chenhui edit end

    tmp<scalarField> intFld = patchInternalField();

    const mappingPotential& nbrField =
    refCast<const mappingPotential>
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            neighbourFieldName_
        )
    );

    // Swap to obtain full local values of neighbour internal field
//***Chenhui edit start
//***Note: OpenFOAM and foam-extend has different formats
//    scalarField nbrIntFld = nbrField.patchInternalField();
//    mapDistribute::distribute
//    (
//        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
//        distMap.schedule(),
//        distMap.constructSize(),
//        distMap.subMap(),           // what to send
//        distMap.constructMap(),     // what to receive
//        nbrIntFld
//    );
     scalarField nbrIntFld = nbrField.patchInternalField();
     mpp.distribute(nbrIntFld);
//***Chenhui edit end
//***Chenhui edit start
//***Note: OpenFOAM and foam-extend has different formats
//    scalarField nbrSigma = nbrField.sigma();
//    mapDistribute::distribute
//    (
//        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
//        distMap.schedule(),
//        distMap.constructSize(),
//        distMap.subMap(),           // what to send
//        distMap.constructMap(),     // what to receive
//        nbrSigma
//    );
      scalarField nbrSigma = nbrField.sigma();
      mpp.distribute(nbrSigma);
//***Chenhui edit end
//***Chenhui edit start
//***Note: OpenFOAM and foam-extend have different formats
//    vectorField nbrEField = nbrField.EField();
//    mapDistribute::distribute
//    (
//        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
//        distMap.schedule(),
//        distMap.constructSize(),
//        distMap.subMap(),           // what to send
//        distMap.constructMap(),     // what to receive
//        nbrEField
//    );
     vectorField nbrEField = nbrField.EField();
     mpp.distribute(nbrEField);
//***Chenhui edit end
//***Chenhui edit start
//***Note: OpenFOAM and foam-extend have different formats
//    scalarField nbrEpsilonDelta = nbrField.epsilon()*nbrPatch.deltaCoeffs();
//    mapDistribute::distribute
//    (
//        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
//        distMap.schedule(),
//        distMap.constructSize(),
//        distMap.subMap(),           // what to send
//        distMap.constructMap(),     // what to receive
//        nbrEpsilonDelta
//    );
     scalarField nbrEpsilonDelta = nbrField.epsilon()*nbrPatch.deltaCoeffs();
     mpp.distribute(nbrEpsilonDelta);
//***Chenhui edit end

    tmp<scalarField> myEpsilonDelta = epsilon()*patch().deltaCoeffs();


    if (neighbourPatchType_ == "dielectric" || neighbourPatchType_ == "otherDielectric"){
        this->refValue() = nbrIntFld;
        this->refGrad() = 0;
        this->valueFraction() = nbrEpsilonDelta / (nbrEpsilonDelta + myEpsilonDelta());  
    }
    else if (neighbourPatchType_ == "plasma"){
        this->refValue() =  nbrIntFld;
        this->refGrad() = 0;
        this->valueFraction() = nbrEpsilonDelta / (nbrEpsilonDelta + myEpsilonDelta());   
    }
    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::mappingPotential::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName")<< neighbourPatchName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchType")<< neighbourPatchType_
        << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		mappingPotential
	);
} // End namespace Foam
// ************************************************************************* //
