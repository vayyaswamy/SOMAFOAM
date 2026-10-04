/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "iterativeCoupledPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::iterativeCoupledPotentialFvPatchScalarField::permittivity
(
    const fvPatch& p
)
{
    const fvMesh& mesh = p.boundaryMesh().mesh();

    const word name
    (
        mesh.foundObject<volScalarField>("epsilonEff")
      ? "epsilonEff"
      : "epsilon"
    );

    return p.patchInternalField
    (
        mesh.lookupObject<volScalarField>(name).internalField()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iterativeCoupledPotentialFvPatchScalarField::
iterativeCoupledPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    remoteFieldName_(iF.name()),
    surfaceCharge_(false),
    fixesValue_(false)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


Foam::iterativeCoupledPotentialFvPatchScalarField::
iterativeCoupledPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    remoteFieldName_(dict.lookupOrDefault<word>("remoteField", iF.name())),
    surfaceCharge_(dict.lookupOrDefault<Switch>("surfaceCharge", false)),
    fixesValue_(false)
{
    if (!isA<regionCoupleFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "iterativeCoupledPotentialFvPatchScalarField::"
            "iterativeCoupledPotentialFvPatchScalarField(...)",
            dict
        )   << "Patch " << p.name() << " is of type " << p.type()
            << ", not " << regionCoupleFvPatch::typeName
            << exit(FatalIOError);
    }

    fixesValue_ = dict.lookupOrDefault<Switch>
    (
        "fixesValue",
        refCast<const regionCoupleFvPatch>(p).master()
    );

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = 0.0;
    }

    valueFraction() = fixesValue_ ? 1.0 : 0.0;
}


Foam::iterativeCoupledPotentialFvPatchScalarField::
iterativeCoupledPotentialFvPatchScalarField
(
    const iterativeCoupledPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    remoteFieldName_(ptf.remoteFieldName_),
    surfaceCharge_(ptf.surfaceCharge_),
    fixesValue_(ptf.fixesValue_)
{}


Foam::iterativeCoupledPotentialFvPatchScalarField::
iterativeCoupledPotentialFvPatchScalarField
(
    const iterativeCoupledPotentialFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    remoteFieldName_(ptf.remoteFieldName_),
    surfaceCharge_(ptf.surfaceCharge_),
    fixesValue_(ptf.fixesValue_)
{}


Foam::iterativeCoupledPotentialFvPatchScalarField::
iterativeCoupledPotentialFvPatchScalarField
(
    const iterativeCoupledPotentialFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    remoteFieldName_(ptf.remoteFieldName_),
    surfaceCharge_(ptf.surfaceCharge_),
    fixesValue_(ptf.fixesValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::iterativeCoupledPotentialFvPatchScalarField&
Foam::iterativeCoupledPotentialFvPatchScalarField::shadow() const
{
    const regionCoupleFvPatch& nbr =
        refCast<const regionCoupleFvPatch>(patch()).shadow();

    return refCast<const iterativeCoupledPotentialFvPatchScalarField>
    (
        nbr.lookupPatchField<volScalarField, scalar>(remoteFieldName_)
    );
}


Foam::tmp<Foam::scalarField>
Foam::iterativeCoupledPotentialFvPatchScalarField::sigma() const
{
    const regionCoupleFvPatch& rcp =
        refCast<const regionCoupleFvPatch>(patch());

    tmp<scalarField> tsigma(new scalarField(size(), 0.0));
    scalarField& s = tsigma();

    if (surfaceCharge_)
    {
        s += patch().lookupPatchField<volScalarField, scalar>("surfC");
    }

    if (shadow().surfaceCharge())
    {
        const scalarField& nbrSurfC =
            rcp.shadow().lookupPatchField<volScalarField, scalar>("surfC");

        s += rcp.interpolate(nbrSurfC);
    }

    return tsigma;
}


Foam::tmp<Foam::scalarField>
Foam::iterativeCoupledPotentialFvPatchScalarField::shadowInterfaceValue() const
{
    const scalarField& nbrValue = shadow();

    return refCast<const regionCoupleFvPatch>(patch()).interpolate(nbrValue);
}


void Foam::iterativeCoupledPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (fixesValue_)
    {
        // The interface potential refValue() is set by the solver
        valueFraction() = 1.0;
    }
    else
    {
        // Gauss's law with the gradient the other side has on the
        // interface: eps dPhi/dn = sigma - epsNbr dPhiNbr/dnNbr
        const regionCoupleFvPatch& rcp =
            refCast<const regionCoupleFvPatch>(patch());

        const fvPatchScalarField& nbrField = shadow();

        const scalarField nbrFlux
        (
            rcp.interpolate
            (
                permittivity(rcp.shadow())*nbrField.snGrad()
            )
        );

        valueFraction() = 0.0;
        refGrad() = (sigma() - nbrFlux)/permittivity(patch());
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::iterativeCoupledPotentialFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("surfaceCharge")
        << surfaceCharge_ << token::END_STATEMENT << nl;
    os.writeKeyword("fixesValue")
        << fixesValue_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        iterativeCoupledPotentialFvPatchScalarField
    );
}

// ************************************************************************* //
