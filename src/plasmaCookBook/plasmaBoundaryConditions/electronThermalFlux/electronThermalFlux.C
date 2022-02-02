/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electronThermalFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronThermalFlux::
electronThermalFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  scalar1Data_(0.0)
{}


Foam::electronThermalFlux::
electronThermalFlux
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
        fvPatchField<vector>::operator=(vector(0.0, 0.0, 0.0));
    }
}


Foam::electronThermalFlux::
electronThermalFlux
(
    const electronThermalFlux& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    scalar1Data_(ptf.scalar1Data_)
{}


Foam::electronThermalFlux::
electronThermalFlux
(
    const electronThermalFlux& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    scalar1Data_(ptf.scalar1Data_)
{}


Foam::electronThermalFlux::
electronThermalFlux
(
    const electronThermalFlux& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    scalar1Data_(ptf.scalar1Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::electronThermalFlux::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::electronThermalFlux::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::electronThermalFlux::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //Info << "electronThermalFlux " << endl;

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

    //const fvPatchField<scalar>& rhof=
    //    patch().lookupPatchField<volScalarField, scalar>("rho");

        vectorField n = patch().nf();

        //Info << "nf = " << n << endl;

        //Info << "Tef = " << Tef << endl;

    	scalarField Enorm = Ef&n;

        vectorField temp = 0.0*n;

        scalarField Finorm = Fif&n;

       // Info << "Finorm " << Finorm << endl;

    forAll(temp, facei)
    {
    	//label faceCelli = patch().faceCells()[facei];

        //Info << "Finorm[facei] " << Finorm[facei] << endl; 
        
		if(Enorm[facei] > 0.0)
		{	
			//Info << "Tef = " << Tef[facei] << endl; 
			//Info << "Nef = " << Nef[facei] << endl;
		    //temp[facei]=(rhof[facei]*(0.25*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei]) 
			//		- scalar1Data_*rhof[facei]*(Fi[faceCelli]/Nef[facei]));
            temp[facei]=((0.25*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei]) 
                    - scalar1Data_*(Finorm[facei]/Nef[facei])*n[facei]);
        }
        else
        {
		    temp[facei]=((0.25*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei]));
        }
	}
    operator == (temp);
    //Info << "temp = " << temp << endl;
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::electronThermalFlux::write(Ostream& os) const
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
        electronThermalFlux
    );
}


// ************************************************************************* //
