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

	//Info << "*this =  " << *this << endl;	
    fvPatchScalarField::evaluate();

    //Info << "*this =  " << *this << endl;
}


void Foam::plasmaDielectricEpsilonFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Info << "*this 1 =  " << *this << endl;
    *this == calcEpsilon(*this, shadowPatchField());
    //Info << "*this 2 =  " << *this << endl;
}


Foam::tmp<Foam::scalarField> 
Foam::plasmaDielectricEpsilonFvPatchScalarField::calcEpsilon
(
    const plasmaDielectricRegionCoupleBase& owner,
    const plasmaDielectricRegionCoupleBase& neighbour
) const
{
    //Info << "calcEpsilon " << endl;
    const fvPatch& p = owner.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);



    //const volScalarField& PhiOwn = this->db().lookupObject<volScalarField>("Phi"); 
    //const volScalarField& PhiOwnOld  = PhiOwn.oldTime();

    //scalar patchIndex = p.index();

    //const scalarField PhiOwnOldPatch = PhiOwnOld.boundaryField()[patchIndex];

    //scalar neighborPatchIndex = neighbour.patch().index();

    //Info << "PhiOwn = " << PhiOwn << endl;

    //const scalarField PhiNeiOldPatch = 

    //const coupledPotentialFvPatchScalarField& PhiwOwn =
    //    dynamic_cast<const coupledPotentialFvPatchScalarField&>
    //    (
    //        PhiOwnOld.boundaryField()[patchIndex]
    //    );

    const coupledPotentialFvPatchScalarField& TwOwn =
        dynamic_cast<const coupledPotentialFvPatchScalarField&>
        (
            p.lookupPatchField<volScalarField, scalar>("Phi")
        );


    //Info << "PhiwOwn = " << TwOwnOld << endl ;

    //Info << "PhiwNei = " << PhiwOwn.shadowPatchField().patchInternalField() << endl;

    //Info << "TwOwn = " << TwOwn << endl;
    //scalarField& k = owner;
    //const scalarField& fOwn = owner.originalPatchField();
    //Info << "here " << endl;
    const scalarField fOwn = neighbour.shadowPatchField().patchInternalField();
    
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Sc(p.size(), 0.0);

    if (TwOwn.surfaceCharge())
    {
        //Info << "I have surface charge " << endl;
        Sc += p.lookupPatchField<volScalarField, scalar>("surfC");
    }

    {
        Field<VectorN<scalar, 4> > lData
        (
            neighbour.size(),
            pTraits<VectorN<scalar, 4> >::zero
        );

        //Info << "here 2 " << endl;
        const scalarField lfNei = owner.shadowPatchField().patchInternalField();
        //Info << "there 2 " << endl;
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

        //Info << "lData = " << lData << endl;
        //Info << "iData = " << iData << endl;

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

    //Info << "1 - p " << (1.0 - p.weights())*mld.magDelta(p.index()) << endl;




    tmp<scalarField> kTmp(new scalarField(p.size()));
    scalarField& k = kTmp();

    //Info << "Sc = " << Sc << endl;

    //k = kOwn*(kNei)/p.deltaCoeffs()/(kOwn + kNei) + 0.5*Sc/p.deltaCoeffs()*(kNei - kOwn)/(kOwn + kNei)/stabilise((TcOwn - TcNei), SMALL);

    //k = kOwn*(kNei)/p.deltaCoeffs()/(kOwn + kNei)*(1.0 + Sc*(kNei/kOwn - 1.0)/kNei/stabilise((TcOwn - TcNei), SMALL));

    k = kOwn*(kNei)/p.deltaCoeffs()/(kOwn + kNei);

    //Info << "TcOwn oldtime = " << p.lookupPatchField<volScalarField, scalar>("Phi").oldTime() << endl;

    //Info << "TwOwn oldtime = " << TwOwn.oldTime() << endl;

    //Info << "Exit calcEpsilon " << endl ;

    // Do interpolation
    harmonic<scalar> interp(mesh);
    const scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = kOwn*(kNei)/p.deltaCoeffs()/(kOwn + kNei);

    //forAll (k, facei)
    //{
    //    k[facei] = max(min(k[facei], 1e2*kHarm[facei]), 1e-2*kHarm[facei]);
    //}

    //Info << "epsilon = " << k << " Sc = " << Sc << " k ratio = " << kNei/kOwn << endl;
    
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
    //Info << "In function calcPotential " << endl;
    const fvPatch& p = TwOwn.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);

    //const scalarField& fOwn = ownerEps.originalPatchField();

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

    //Info << "kOwn = " << kOwn << endl;

    //Info << "fOwn inside calcPotential = " << fOwn << endl;

    //Info << "fNei inside calcPotential = " << fNei << endl;



    //Info << "magDelta = " << mld.magDelta(p.index()) << endl;

    //Info << "weights = " << p.weights() << endl;

    //Info << "delta = " << p.deltaCoeffs() << endl;

    //Info << "kNei = " << kNei << endl;

    //Info << "TcOwn = " << TcOwn << endl;

    //Info << "TcNei = " << TcNei << endl;

    //scalarField Tw1 = (-Sc/fOwn*deltaOwn + TcOwn + fNei/fOwn*TcNei*deltaOwn/deltaNei)/(1.0 + fNei/fOwn*deltaOwn/deltaNei);

    //Info << "Tw1 = " << Tw1 << endl;

    //scalarField a = pos((Tw1 - TcOwn)*(Tw1 - TcNei));

    //Info << "a = " << a << endl ;

    //scalarField Tw2 = (Sc/fOwn*deltaOwn + TcOwn + fNei/fOwn*TcNei*deltaOwn/deltaNei)/(1.0 + fNei/fOwn*deltaOwn/deltaNei);

    //Info << "Tw2  = " << Tw2 << endl;

    //Tw = a*Tw2 + (1-a)*Tw1;

    //Tw = Tw1 ;

    //Info << "Sc = " << Sc << endl;

    //Tw = (kOwn*TcOwn + kNei*TcNei + Sc)/(kOwn + kNei);

    Tw = kOwn*TcOwn*(1.0 + kNei/kOwn*TcNei/TcOwn + Sc/kOwn/TcOwn)/(kOwn + kNei);

    //Info << "Tw = " << Tw << endl;

    //scalarField q1 = -(TcOwn - Tw)*kOwn;

    //scalarField q2 = -(TcNei - Tw)*kNei;

    //scalarField q3 = (TcNei - TcOwn)*ownerEps*p.deltaCoeffs();

    //Info << "q1 = " << q1 << endl;

    //Info << "q2 = " << q2 << endl;

    //Info << "q3 = " << q3 << endl;

    //Info << "q2 + q1 " << q2 + q1 << endl;

    //Info << "E1 = " << (Tw-TcOwn)/(1.0 - p.weights())/mld.magDelta(p.index()) << endl;

    //Info << "E2 = " << (TcNei-Tw)/(p.weights())/mld.magDelta(p.index()) << endl;

    //Info << "average " << 0.5*(q1 + q2) << endl;

    //Info << "ownerEps = " << ownerEps << endl;

    return TwTmp;


    //const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    //const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    //TwOwn ==
    //    (kOwn*TcOwn + kNei*TcNei + Sc)
    //   /((kOwn + kNei));

    //TwOwn ==
    //    (kOwn*TcOwn + kNei*TcNei)
    //   /((kOwn + kNei));

    //Info << "E = " << kOwn*(TcOwn - TwOwn)/(1.0 - p.weights())/mld.magDelta(p.index()) << endl;
    //Info << "Sc = " << Sc << endl;

    //Info << "calcPotential end" << endl;

    //TwOwn.fvPatchScalarField::updateCoeffs();
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
