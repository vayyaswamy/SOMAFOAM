/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "eEnergyHGLFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eEnergyHGLFlux::
eEnergyHGLFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  scalar1Data_(0.0),
  scalar2Data_(0.0),
  scalar3Data_(0.0)
{}


Foam::eEnergyHGLFlux::
eEnergyHGLFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    scalar1Data_(readScalar(dict.lookup("sec"))),
    scalar2Data_(readScalar(dict.lookup("refCoeff"))),
    scalar3Data_(readScalar(dict.lookup("secTemp")))
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


Foam::eEnergyHGLFlux::
eEnergyHGLFlux
(
    const eEnergyHGLFlux& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_),
    scalar3Data_(ptf.scalar3Data_)
{}


Foam::eEnergyHGLFlux::
eEnergyHGLFlux
(
    const eEnergyHGLFlux& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_),
    scalar3Data_(ptf.scalar3Data_)
{}


Foam::eEnergyHGLFlux::
eEnergyHGLFlux
(
    const eEnergyHGLFlux& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_),
    scalar3Data_(ptf.scalar3Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::eEnergyHGLFlux::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::eEnergyHGLFlux::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::eEnergyHGLFlux::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volVectorField& Fi =
        db().objectRegistry::lookupObject<volVectorField>("ionFlux");

    const fvPatchField<scalar>& Nef =
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te");

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<scalar>& muef=
        patch().lookupPatchField<volScalarField, scalar>("mu_electron");

        vectorField n = patch().nf();

        label patchi = this->patch().index();
        const fvPatchVectorField& FiT = Fi.boundaryField()[patchi];

    	scalarField Enorm = Ef&n;
		scalarField Nf = scalar1Data_*(FiT&n)/muef/(Enorm+0.001);
        vectorField temp = 0.0*n;        

    forAll(temp, facei)
    {
    	label faceCelli = patch().faceCells()[facei];

		if(Enorm[facei] > 0.0)
		{	
		    temp[facei]=((((1-scalar2Data_)/(1+scalar2Data_))*3.4516212999999994e-23*((Tef[facei]*Nef[facei])-(Nf[facei]*scalar3Data_))*((0.50*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei])+muef[facei]*Ef[facei])) 
					- scalar1Data_*scalar3Data_*3.4516212999999994e-23*Fi[faceCelli]);
        }
        else
        {        
		    temp[facei]=(((1-scalar2Data_)/(1+scalar2Data_))*3.4516212999999994e-23*(Tef[facei]*Nef[facei])*((0.50*Foam::sqrt(3.8595300515e7*Tef[facei])*n[facei])+muef[facei]*Ef[facei]));
        }
	}

    operator == (temp);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::eEnergyHGLFlux::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("sec")
        << scalar1Data_ << token::END_STATEMENT << nl;
    os.writeKeyword("refCoeff")
        << scalar2Data_ << token::END_STATEMENT << nl;
    os.writeKeyword("secTemp")
        << scalar3Data_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        eEnergyHGLFlux
    );
}


// ************************************************************************* //
