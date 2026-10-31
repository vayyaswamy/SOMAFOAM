/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "ionVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionVelocity::
ionVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
    fieldName_(iF.name()),
  scalar1Data_(0.0)
{}


Foam::ionVelocity::
ionVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    fieldName_(iF.name()),
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


Foam::ionVelocity::
ionVelocity
(
    const ionVelocity& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    scalar1Data_(ptf.scalar1Data_)
{}


Foam::ionVelocity::
ionVelocity
(
    const ionVelocity& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    fieldName_(ptf.fieldName_),
    scalar1Data_(ptf.scalar1Data_)
{}


Foam::ionVelocity::
ionVelocity
(
    const ionVelocity& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    scalar1Data_(ptf.scalar1Data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::ionVelocity::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::ionVelocity::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::ionVelocity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    Info << "it came till here" << endl;
    const volScalarField& Nef =
        db().objectRegistry::lookupObject<volScalarField>("Arp1");
    const volVectorField& vef =
        db().objectRegistry::lookupObject<volVectorField>(fieldName_);
    const fvPatchField<scalar>& Nef_bound =
        patch().lookupPatchField<volScalarField, scalar>("Arp1");
    const fvPatchField<vector>& vef_bound=
        patch().lookupPatchField<volVectorField, scalar>(fieldName_);

    vectorField n = patch().nf();

    label patchi = this->patch().index();

    vectorField temp = 0.0*n;
    

    forAll(temp, facei)
    {
        label faceCelli = patch().faceCells()[facei];
        int a= pos(vef[facei]&n[facei]);
		temp[facei]=a*vef[facei];//+(0.25*Foam::sqrt(8*1.38e-23*etemperat/(22/7)/9.1e-31)*n[facei]) ;
	}

    operator == (temp);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::ionVelocity::write(Ostream& os) const
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
        ionVelocity
    );
}


// ************************************************************************* //