/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "eEnergyThermalMFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eEnergyThermalMFlux::
eEnergyThermalMFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  scalar1Data_(0.0),
  scalar2Data_(0.0)
  {}


Foam::eEnergyThermalMFlux::
eEnergyThermalMFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    scalar1Data_(readScalar(dict.lookup("sec"))),
    scalar2Data_(readScalar(dict.lookup("secTemp")))
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


Foam::eEnergyThermalMFlux::
eEnergyThermalMFlux
(
    const eEnergyThermalMFlux& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_)
{}


Foam::eEnergyThermalMFlux::
eEnergyThermalMFlux
(
    const eEnergyThermalMFlux& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_)
{}


Foam::eEnergyThermalMFlux::
eEnergyThermalMFlux
(
    const eEnergyThermalMFlux& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::eEnergyThermalMFlux::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::eEnergyThermalMFlux::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::eEnergyThermalMFlux::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volVectorField& Fi =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux");

    const fvPatchField<vector>& Jef =
        patch().lookupPatchField<volVectorField, vector>("J_electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te");

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

        vectorField n = patch().nf();

    	scalarField Enorm = Ef&n;
        vectorField temp = 0.0*n;        

    forAll(temp, facei)
    {
    	label faceCelli = patch().faceCells()[facei];

		if(Enorm[facei] > 0.0)
		{	
		    temp[facei]=(2.76129704e-23*(Tef[facei]*Jef[facei])
					- scalar1Data_*scalar2Data_*2.76129704e-23*Fi[faceCelli]);
        }
        else
        {        
		    temp[facei]=(2.76129704e-23*(Tef[facei]*Jef[facei]));
        }
	}

    operator == (temp);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::eEnergyThermalMFlux::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("sec")
        << scalar1Data_ << token::END_STATEMENT << nl;
    os.writeKeyword("secTemp")
        << scalar2Data_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        eEnergyThermalMFlux
    );
}


// ************************************************************************* //
