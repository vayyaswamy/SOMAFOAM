/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "ionHGLFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"


#include "linear.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionHGLFlux::
ionHGLFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
	fixedValueFvPatchVectorField(p, iF),
	mobName_("mobility"),
	scalar1Data_(0.0),
	scalar2Data_(0.0),
	scalar3Data_(0.0)    
{}


Foam::ionHGLFlux::
ionHGLFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
	mobName_(dict.lookupOrDefault<word>("mobility","mu_Arp1")),
    scalar1Data_(readScalar(dict.lookup("refCoeff"))),
    scalar2Data_(readScalar(dict.lookup("charge"))),
    scalar3Data_(readScalar(dict.lookup("molWeight")))
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


Foam::ionHGLFlux::
ionHGLFlux
(
    const ionHGLFlux& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    mobName_(ptf.mobName_),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_),
    scalar3Data_(ptf.scalar3Data_)
{}


Foam::ionHGLFlux::
ionHGLFlux
(
    const ionHGLFlux& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    mobName_(ptf.mobName_),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_),
    scalar3Data_(ptf.scalar3Data_)
{}


Foam::ionHGLFlux::
ionHGLFlux
(
    const ionHGLFlux& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    mobName_(ptf.mobName_),
    scalar1Data_(ptf.scalar1Data_),
    scalar2Data_(ptf.scalar2Data_),
    scalar3Data_(ptf.scalar3Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::ionHGLFlux::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::ionHGLFlux::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::ionHGLFlux::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatchField<scalar>& Tif=
        patch().lookupPatchField<volScalarField, scalar>("Tion");

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<scalar>& muif=
        patch().lookupPatchField<volScalarField, scalar>(mobName_);

    //const fvPatchField<scalar>& rhof=
    //    patch().lookupPatchField<volScalarField, scalar>("rho");

        vectorField n = patch().nf();

    	scalarField Enorm = scalar2Data_*(Ef&n);
        vectorField temp = 0.0*n;        

    forAll(temp, facei)
    {
		if(Enorm[facei] > 0.0)
		{	
		    temp[facei]=(((1-scalar1Data_)/(1+scalar1Data_))*((0.50*Foam::sqrt(21172.59783*Tif[facei]/scalar3Data_)*n[facei])+muif[facei]*Ef[facei]));
        }
        else
        {        
		    temp[facei]=(((1-scalar1Data_)/(1+scalar1Data_))*((0.50*Foam::sqrt(21172.59783*Tif[facei]/scalar3Data_)*n[facei])-muif[facei]*Ef[facei]));
        }
	}

    operator == (temp);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::ionHGLFlux::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("mobility")
        << mobName_ << token::END_STATEMENT << nl;    
    os.writeKeyword("refCoeff")
        << scalar1Data_ << token::END_STATEMENT << nl;
    os.writeKeyword("charge")
        << scalar2Data_ << token::END_STATEMENT << nl;
    os.writeKeyword("molWeight")
        << scalar3Data_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        ionHGLFlux
    );
}


// ************************************************************************* //
