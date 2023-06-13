/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaExtendedPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

namespace Foam
{
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::scalarField
Foam::plasmaExtendedPotential::voltage() const
{
    // Extract the dictionary from the database
    IOdictionary& voltageDict = const_cast<IOdictionary&>(db().lookupObject<IOdictionary>
    (
        "voltageDict"
    ));

    // Extracting scalar value
    scalar voltage(voltageDict.lookupOrDefault<scalar>("voltage",0.0));

    return tmp<scalarField>
    (
        new scalarField(this->size(), voltage)
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::plasmaExtendedPotential::
plasmaExtendedPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::plasmaExtendedPotential::
plasmaExtendedPotential
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}

Foam::plasmaExtendedPotential::
plasmaExtendedPotential
(
    const plasmaExtendedPotential& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::plasmaExtendedPotential::
plasmaExtendedPotential
(
    const plasmaExtendedPotential& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::plasmaExtendedPotential::
plasmaExtendedPotential
(
    const plasmaExtendedPotential& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::plasmaExtendedPotential::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    plasmaExtendedPotential::autoMap(m);
}


void Foam::plasmaExtendedPotential::
rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    plasmaExtendedPotential::rmap(ptf, addr);
}

void Foam::plasmaExtendedPotential::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(voltage());

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::plasmaExtendedPotential::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, plasmaExtendedPotential);

}
// ************************************************************************* //
