/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electronThermalVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronThermalVelocity::
electronThermalVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  scalar1Data_(0.0)
{}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    scalar1Data_(readScalar(dict.lookup("sec")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }
}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const electronThermalVelocity& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    scalar1Data_(ptf.scalar1Data_)
{}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const electronThermalVelocity& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    scalar1Data_(ptf.scalar1Data_)
{}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const electronThermalVelocity& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    scalar1Data_(ptf.scalar1Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::electronThermalVelocity::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::electronThermalVelocity::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::electronThermalVelocity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    Info << "it came till here" << endl;
    const volVectorField& Fi =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux");

    const volScalarField& TefCell =
        db().objectRegistry::lookupObject<volScalarField>("eTemp");

    const fvPatchField<scalar>& Nef =
        patch().lookupPatchField<volScalarField, scalar>("electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("eTemp");

    vectorField n = patch().nf();

    label patchi = this->patch().index();

    vectorField temp = 0.0*n;
    

    forAll(temp, facei)
    {
        label faceCelli = patch().faceCells()[facei];
		temp[facei]=(0.25*Foam::sqrt(8*1.38e-23*Tef[facei]/(22/7)/9.1e-31)*n[facei]) ;
	}

    operator == (temp);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::electronThermalVelocity::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("sec")
        << scalar1Data_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        electronThermalVelocity
    );
}


// ************************************************************************* //