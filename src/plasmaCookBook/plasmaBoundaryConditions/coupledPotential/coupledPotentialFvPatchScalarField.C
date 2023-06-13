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
#include "magLongDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledPotentialFvPatchScalarField::coupledPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF),
	surfaceCharge_(false)
{}


coupledPotentialFvPatchScalarField::coupledPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict),
    surfaceCharge_(readBool(dict.lookup("surfaceCharge")))
{}


coupledPotentialFvPatchScalarField::coupledPotentialFvPatchScalarField
(
    const coupledPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper),
    surfaceCharge_(ptf.surfaceCharge_)
{}


coupledPotentialFvPatchScalarField::coupledPotentialFvPatchScalarField
(
    const coupledPotentialFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF),
    surfaceCharge_(ptf.surfaceCharge_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Return a named shadow patch field
const coupledPotentialFvPatchScalarField&
coupledPotentialFvPatchScalarField::shadowPatchField() const
{
    return dynamic_cast
    <
        const coupledPotentialFvPatchScalarField&
    >
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}

// Return neighbour field given internal cell data
tmp<scalarField>
coupledPotentialFvPatchScalarField::patchNeighbourField() const
{
    return regionCouplingFvPatchScalarField::patchNeighbourField
    (
        remoteFieldName()
    );
}


tmp<scalarField>
coupledPotentialFvPatchScalarField::Phiw() const
{
    return *this;
}


tmp<scalarField>
coupledPotentialFvPatchScalarField::Phic() const
{
    return patchInternalField();
}


void coupledPotentialFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{

    const fvPatchScalarField& epspf =
        lookupPatchField<volScalarField, scalar>("epsilon");

    const plasmaDielectricRegionCoupleBase& eps =
        dynamic_cast<const plasmaDielectricRegionCoupleBase&>(epspf);

    //Info << "starting calcPotential " << endl;



    *this == eps.calcPotential(*this, shadowPatchField(), eps);

    //Info << "initEvaluate " << endl;
}


void coupledPotentialFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    //Info << "about to evaluate " << endl;
    //Info << "Shadow = " << shadowPatchField() << endl; 
    fvPatchScalarField::evaluate();
    //Info << "evaluate " << endl;
}


void coupledPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    //Info << "About to do updateCoeffs " << endl;

    const plasmaDielectricRegionCoupleBase& eps = 
        refCast<const plasmaDielectricRegionCoupleBase>
        (
            lookupPatchField<volScalarField, scalar>("epsilon")
        );

    //Info << "About to set this " << endl;

    *this == eps.calcPotential(*this, shadowPatchField(), eps);

    fvPatchScalarField::updateCoeffs();
    //Info << "Done with updateCoeffs " << endl;
}


tmp<scalarField> coupledPotentialFvPatchScalarField::source() const
{
    const fvPatch& p = patch();

    const scalarField PhicOwn = Phic();
    const scalarField PhicNei =
        regionCouplePatch().interpolate(shadowPatchField().Phic());
    const scalarField Phiw = this->Phiw();

    const plasmaDielectricRegionCoupleBase& Eps =
        dynamic_cast<const plasmaDielectricRegionCoupleBase&>
        (
            p.lookupPatchField<volScalarField, scalar>("epsilon")
        );

    const scalarField eps = Eps.kw()*p.deltaCoeffs();
    const scalarField epsOwn = Eps.kc();

    //Info << "epsOwn = " << epsOwn << endl;

    //Info << "Eps.kw() = " << Eps.kw() << endl;

    //Info << "eps =  " << Eps.kw()*p.deltaCoeffs() << endl;

    //Info << "source = " << epsOwn*(PhicOwn - Phiw) + eps*(PhicNei - PhicOwn) << endl;

    return epsOwn*(PhicOwn - Phiw) + eps*(PhicNei - PhicOwn);

}


void coupledPotentialFvPatchScalarField::manipulateMatrix
(
    fvScalarMatrix& matrix
)
{
    const fvPatch& p = patch();
    const scalarField& magSf = p.magSf();
    const labelList& cellLabels = p.faceCells();
    scalarField& source = matrix.source();

    //Info << "About to call source " << endl;

    scalarField s = this->source();

    //Info << "s = " << s << endl;

    forAll(cellLabels, i)
    {
        source[cellLabels[i]] += s[i]*magSf[i];
    }
    //Info << "Done with manipulateMatrix" << endl;
}


void coupledPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("surfaceCharge") << surfaceCharge_ << token::END_STATEMENT << nl;
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    coupledPotentialFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
