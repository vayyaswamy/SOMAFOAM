/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "boundaryConditions/mappingPotential/mappingPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"


namespace Foam
{
mappingPotential::mappingPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nbrFieldName_("undefined-nbrFieldName"),
    nbrPatchName_("undefined-nbrPatchName")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

mappingPotential::mappingPotential
(
    const mappingPotential& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    nbrFieldName_(psf.nbrFieldName_),
    nbrPatchName_(psf.nbrPatchName_)
{}

mappingPotential::mappingPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nbrFieldName_(dict.lookup("nbrFieldName")),
    nbrPatchName_(dict.lookup("nbrPatchName"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
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

mappingPotential::mappingPotential
(
    const mappingPotential& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    nbrFieldName_(psf.nbrFieldName_),
    nbrPatchName_(psf.nbrPatchName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const fvPatchScalarField&
mappingPotential::sigma() const
{
    const fvPatchField<scalar>& tsurfC =
        patch().lookupPatchField<volScalarField, scalar>("surfC");
    return tsurfC;
}

const fvPatchVectorField&
mappingPotential::EField() const
{
    const fvPatchField<vector>& tE =
        patch().lookupPatchField<volVectorField, vector>("EField");
    return tE;
}

const fvPatchScalarField&
mappingPotential::epsilon() const
{
    const fvPatchField<scalar>& tEps =
        patch().lookupPatchField<volScalarField, scalar>("epsilon");
    return tEps;
}

const fvPatchScalarField&
mappingPotential::phi() const
{
    const fvPatchField<scalar>& tphi =
        patch().lookupPatchField<volScalarField, scalar>("EPotential");
    return tphi;
}
void mappingPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the directMappedPatchBase
    const directMappedPatchBase& mpp =
    refCast<const directMappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
    refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];
    const mappingPotential&
        nbrField = refCast
            <const mappingPotential>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(nbrFieldName_)
            );
     scalarField nbrIntFld = nbrField.patchInternalField();
     mpp.distribute(nbrIntFld);

     scalarField nbrSigma = nbrField.sigma();
     mpp.distribute(nbrSigma);

     scalarField nbrEpsilonDelta = nbrField.epsilon()*nbrPatch.deltaCoeffs();
     mpp.distribute(nbrEpsilonDelta);

     vectorField nbrEField = nbrField.EField();
     mpp.distribute(nbrEField);
    tmp<scalarField> myEpsilonDelta = epsilon()*patch().deltaCoeffs();

    this->refValue() =  nbrIntFld;
    this->refGrad() = 0;
    this->valueFraction() = nbrEpsilonDelta / (nbrEpsilonDelta + myEpsilonDelta());   
    mixedFvPatchScalarField::updateCoeffs();

    UPstream::msgType() = oldTag;
}

void mappingPotential::write
(
    Ostream& os
) const
{

    mixedFvPatchScalarField::write(os);
      os.writeKeyword("nbrFieldName") << nbrFieldName_ << token::END_STATEMENT << nl;
      os.writeKeyword("nbrPatchName") << nbrPatchName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
makePatchTypeField
(
    fvPatchScalarField,
    mappingPotential
);

} 
// End namespace Foam
// ************************************************************************* //
