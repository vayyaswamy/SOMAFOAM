/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "plasmaDielectricEpsilonFvPatchScalarField.H"
#include "plasmaDielectricEpsilonSlaveFvPatchScalarField.H"
#include "coupledPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "harmonic.H"
#include "VectorN.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaDielectricEpsilonFvPatchScalarField::
plasmaDielectricEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    plasmaDielectricRegionCoupleBase(p, iF)
{}


Foam::plasmaDielectricEpsilonFvPatchScalarField::
plasmaDielectricEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    plasmaDielectricRegionCoupleBase(p, iF, dict)
{}


Foam::plasmaDielectricEpsilonFvPatchScalarField::
plasmaDielectricEpsilonFvPatchScalarField
(
    const plasmaDielectricEpsilonFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    plasmaDielectricRegionCoupleBase(ptf, p, iF, mapper)
{}


Foam::plasmaDielectricEpsilonFvPatchScalarField::
plasmaDielectricEpsilonFvPatchScalarField
(
    const plasmaDielectricEpsilonFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    plasmaDielectricRegionCoupleBase(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaDielectricEpsilonFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{

    fvPatchScalarField::evaluate();

}


void Foam::plasmaDielectricEpsilonFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    *this == calcEpsilon(*this, shadowPatchField());
}


Foam::tmp<Foam::scalarField>
Foam::plasmaDielectricEpsilonFvPatchScalarField::calcEpsilon
(
    const plasmaDielectricRegionCoupleBase& owner,
    const plasmaDielectricRegionCoupleBase& neighbour
) const
{
    const fvPatch& p = owner.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);


    const coupledPotentialFvPatchScalarField& TwOwn =
        dynamic_cast<const coupledPotentialFvPatchScalarField&>
        (
            p.lookupPatchField<volScalarField, scalar>("Phi")
        );


    const scalarField fOwn = neighbour.shadowPatchField().patchInternalField();

    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Sc(p.size(), 0.0);

    if (TwOwn.surfaceCharge())
    {
        Sc += p.lookupPatchField<volScalarField, scalar>("surfC");
    }

    {
        Field<VectorN<scalar, 4> > lData
        (
            neighbour.size(),
            pTraits<VectorN<scalar, 4> >::zero
        );

        const scalarField lfNei = owner.shadowPatchField().patchInternalField();
        scalarField lTcNei = TwOwn.shadowPatchField().patchInternalField();

        forAll (lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
        }

        if (TwOwn.shadowPatchField().surfaceCharge())
        {
            const scalarField& lQrNei =
                owner.lookupShadowPatchField<volScalarField, scalar>("surfC");
            const scalarField& lTwNei = TwOwn.shadowPatchField();

            forAll (lData, facei)
            {
                lData[facei][2] = lTwNei[facei];
                lData[facei][3] = lQrNei[facei];
            }
        }

        const Field<VectorN<scalar, 4> > iData =
            owner.regionCouplePatch().interpolate(lData);


        forAll (iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
        }

        if (TwOwn.shadowPatchField().surfaceCharge())
        {
            forAll (iData, facei)
            {
                Sc[facei] += iData[facei][3];
            }
        }
    }


    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    const scalarField deltaOwn = (1.0 - p.weights())*mld.magDelta(p.index());

    const scalarField deltaNei = p.weights()*mld.magDelta(p.index());


    tmp<scalarField> kTmp(new scalarField(p.size()));
    scalarField& k = kTmp();


    k = kOwn*(kNei)/p.deltaCoeffs()/(kOwn + kNei);


    // Do interpolation
    harmonic<scalar> interp(mesh);
    const scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = kOwn*(kNei)/p.deltaCoeffs()/(kOwn + kNei);


    return kTmp;

}


Foam::tmp<Foam::scalarField>
Foam::plasmaDielectricEpsilonFvPatchScalarField::calcPotential
(
    const coupledPotentialFvPatchScalarField& TwOwn,
    const coupledPotentialFvPatchScalarField& neighbour,
    const plasmaDielectricRegionCoupleBase& ownerEps
) const
{
    const fvPatch& p = TwOwn.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);


    const scalarField fOwn = ownerEps.patchInternalField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Sc(p.size(), 0.0);

    if (TwOwn.surfaceCharge())
    {
        Sc += p.lookupPatchField<volScalarField, scalar>("surfC");
    }

    {
        Field<VectorN<scalar, 4> > lData
        (
            neighbour.size(),
            pTraits<VectorN<scalar, 4> >::zero
        );

        const scalarField lfNei =
            ownerEps.shadowPatchField().patchInternalField();
        scalarField lTcNei =
            TwOwn.shadowPatchField().patchInternalField();

        forAll (lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
        }

        if (TwOwn.shadowPatchField().surfaceCharge())
        {
            const scalarField& lTwNei = TwOwn.shadowPatchField();
            const scalarField& lQrNei =
                TwOwn.lookupShadowPatchField<volScalarField, scalar>("surfC");

            forAll (lData, facei)
            {
                lData[facei][2] = lTwNei[facei];
                lData[facei][3] = lQrNei[facei];
            }
        }

        const Field<VectorN<scalar, 4> > iData =
            TwOwn.regionCouplePatch().interpolate(lData);

        forAll (iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
        }

        if (TwOwn.shadowPatchField().surfaceCharge())
        {
            forAll (iData, facei)
            {
                Sc[facei] += iData[facei][3];
            }
        }
    }


    harmonic<scalar> interp(mesh);
    scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    const scalarField deltaOwn = (1.0 - p.weights())*mld.magDelta(p.index());

    const scalarField deltaNei = p.weights()*mld.magDelta(p.index());

    tmp<scalarField> TwTmp(new scalarField(TwOwn.Phiw()));
    scalarField& Tw = TwTmp();


    Tw = kOwn*TcOwn*(1.0 + kNei/kOwn*TcNei/TcOwn + Sc/kOwn/TcOwn)/(kOwn + kNei);


    return TwTmp;


}


void Foam::plasmaDielectricEpsilonFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    makePatchTypeField
    (
        fvPatchScalarField,
        plasmaDielectricEpsilonFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
