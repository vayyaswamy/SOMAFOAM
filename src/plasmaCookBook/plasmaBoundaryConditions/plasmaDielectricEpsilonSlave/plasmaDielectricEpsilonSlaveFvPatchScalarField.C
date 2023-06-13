/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaDielectricEpsilonSlaveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::
plasmaDielectricEpsilonSlaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    plasmaDielectricRegionCoupleBase(p, iF)
{}


Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::
plasmaDielectricEpsilonSlaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    plasmaDielectricRegionCoupleBase(p, iF, dict)
{}


Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::
plasmaDielectricEpsilonSlaveFvPatchScalarField
(
    const plasmaDielectricEpsilonSlaveFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    plasmaDielectricRegionCoupleBase(ptf, iF)
{}


Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::
plasmaDielectricEpsilonSlaveFvPatchScalarField
(
    const plasmaDielectricEpsilonSlaveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    plasmaDielectricRegionCoupleBase(ptf, p, iF, mapper)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::plasmaDielectricRegionCoupleBase&
Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::shadowPatchField() const
{
    return dynamic_cast<const plasmaDielectricRegionCoupleBase&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


void Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::plasmaDielectricEpsilonSlaveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //    shadowPatchField().calcEpsilon(*this, shadowPatchField());

    scalarField& k = *this;

    k = calcEpsilon(*this,shadowPatchField());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    makePatchTypeField
    (
        fvPatchScalarField,
        plasmaDielectricEpsilonSlaveFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
