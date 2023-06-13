/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "coupledPotentialFvPatchScalarField.H"
#include "plasmaDielectricRegionCoupleBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaDielectricRegionCoupleBase::plasmaDielectricRegionCoupleBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF)
{}


plasmaDielectricRegionCoupleBase::plasmaDielectricRegionCoupleBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict)
{}


plasmaDielectricRegionCoupleBase::plasmaDielectricRegionCoupleBase
(
    const plasmaDielectricRegionCoupleBase& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper)
{}


plasmaDielectricRegionCoupleBase::plasmaDielectricRegionCoupleBase
(
    const plasmaDielectricRegionCoupleBase& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a named shadow patch field
const plasmaDielectricRegionCoupleBase&
plasmaDielectricRegionCoupleBase::shadowPatchField() const
{
    return dynamic_cast<const plasmaDielectricRegionCoupleBase&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}

tmp<scalarField> plasmaDielectricRegionCoupleBase::forig() const
{
   
    return originalPatchField();
}

tmp<scalarField> plasmaDielectricRegionCoupleBase::kc() const
{
    const fvPatch& p = patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    //Info << "kc func " << patchInternalField() << endl;
   
    return patchInternalField()/(1 - p.weights())/mld.magDelta(p.index());
}


// Return a named shadow patch field
tmp<scalarField> plasmaDielectricRegionCoupleBase::korig() const
{
    const fvPatch& p = patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    return forig()/(1 - p.weights())/mld.magDelta(p.index());
}


// Return a named shadow patch field
tmp<scalarField> plasmaDielectricRegionCoupleBase::kw() const
{
    
    return *this;
    
}


void plasmaDielectricRegionCoupleBase::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();
}


tmp<scalarField>
plasmaDielectricRegionCoupleBase::calcEpsilon
(
    const plasmaDielectricRegionCoupleBase& owner,
    const plasmaDielectricRegionCoupleBase& neighbour
) const
{
    return shadowPatchField().calcEpsilon(owner, neighbour);
}


tmp<scalarField>
plasmaDielectricRegionCoupleBase::calcPotential
(
    const coupledPotentialFvPatchScalarField& owner,
    const coupledPotentialFvPatchScalarField& neighbour,
    const plasmaDielectricRegionCoupleBase& ownerEps
) const
{
    //Info << "owner = " << owner << endl;
    //Info << "neighbour = " << neighbour << endl;
    //Info << "ownerEps = " << ownerEps << endl;
    //Info << "Shadow = " << shadowPatchField() << endl;
    return shadowPatchField().calcPotential(owner, neighbour, ownerEps);
    //Info << "Done here " << endl;
}


void
plasmaDielectricRegionCoupleBase::initInterfaceMatrixUpdate
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes
) const
{
    FatalErrorIn
    (
        "plasmaDielectricRegionCoupleBase::initInterfaceMatrixUpdate"
    )   << abort(FatalError);
}


void
plasmaDielectricRegionCoupleBase::updateInterfaceMatrix
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes
) const
{
    FatalErrorIn
    (
        "plasmaDielectricRegionCoupleBase::updateInterfaceMatrix"
    )   << abort(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    plasmaDielectricRegionCoupleBase
);

} // End namespace Foam

// ************************************************************************* //
