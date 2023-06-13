/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "eEnergyThermalFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eEnergyThermalFlux::
eEnergyThermalFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  scalar1Data_(0.0),
  scalar2Data_(0.0)
  {}


Foam::eEnergyThermalFlux::
eEnergyThermalFlux
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


Foam::eEnergyThermalFlux::
eEnergyThermalFlux
(
    const eEnergyThermalFlux& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_)
{}


Foam::eEnergyThermalFlux::
eEnergyThermalFlux
(
    const eEnergyThermalFlux& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_)
{}


Foam::eEnergyThermalFlux::
eEnergyThermalFlux
(
    const eEnergyThermalFlux& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::eEnergyThermalFlux::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::eEnergyThermalFlux::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::eEnergyThermalFlux::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Info << "eEnergyThermalFlux " << endl;

    //const volVectorField& Fi =
    //    db().objectRegistry::lookupObject<volVectorField>("ionFlux");

    const fvPatchField<scalar>& Nef =
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te");

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

        vectorField n = patch().nf();

    	scalarField Enorm = Ef&n;
        vectorField temp = 0.0*n;  
        scalarField Finorm = Fif&n;      

        //Info << "Finorm " << Finorm << endl;

    forAll(temp, facei)
    {
    	//label faceCelli = patch().faceCells()[facei];

         //Info << "Finorm = " << Finorm[facei] << endl;

		if(Enorm[facei] > 0.0)
		{	

		    temp[facei]=((3.4516212999999994e-23*(Tef[facei]*Nef[facei])*(0.25*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei])) 
					- scalar1Data_*scalar2Data_*3.4516212999999994e-23*Finorm[facei]*n[facei]);
        }
        else
        {        
		    temp[facei]=(3.4516212999999994e-23*(Tef[facei]*Nef[facei])*(0.25*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei]));
        }
	}

    operator == (temp);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::eEnergyThermalFlux::write(Ostream& os) const
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
        eEnergyThermalFlux
    );
}


// ************************************************************************* //
