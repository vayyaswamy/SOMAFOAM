/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electronThermalVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronThermalVelocity::
electronThermalVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  seec_(0),
  Edepend_(true),
  FE_(false),
  beta_(1.0),
  wf_(1.0)
{
    fvPatchField<vector>::operator=(this->patchInternalField());
}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    seec_(readScalar(dict.lookup("seec"))),
    Edepend_(readBool(dict.lookup("Edepend"))),
    FE_(readBool(dict.lookup("field_emission"))),
    beta_(readScalar(dict.lookup("field_enhancement_factor"))),
    wf_(readScalar(dict.lookup("work_function")))
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


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const electronThermalVelocity& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const electronThermalVelocity& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


Foam::electronThermalVelocity::
electronThermalVelocity
(
    const electronThermalVelocity& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    seec_(ptf.seec_),
    Edepend_(ptf.Edepend_),
    FE_(ptf.FE_),
    beta_(ptf.beta_),
    wf_(ptf.wf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::electronThermalVelocity::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::electronThermalVelocity::
rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void Foam::electronThermalVelocity::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField n = patch().nf();

    const fvPatchField<vector>& Fif=
        patch().lookupPatchField<volVectorField, vector>("ionFlux");

    const fvPatchField<scalar>& Tef=
        patch().lookupPatchField<volScalarField, scalar>("Te"); 

    const fvPatchField<vector>& Ef=
        patch().lookupPatchField<volVectorField, vector>("E");

    const fvPatchField<scalar>& Nef =
        patch().lookupPatchField<volScalarField, scalar>("N_electron");

    scalarField Enorm = Ef&n ;

    scalarField Fifnorm = Fif&n;

    scalarField a = pos(Enorm);   

    scalarField b = pos(Fifnorm);

    const scalarField Gamma_se = seec_*(b*Fifnorm);

    scalarField c = pos(Enorm); // to ensure field emission happens only if E-field is pointing inward even if Edepend is set to false

    scalarField Gamma_FE = c*0.0;

    //Info << "FE_ = " << FE_ << endl;

    if (FE_)
    {
        //Info << "Enorm = " << Enorm << endl;

        scalarField vofy = 0.95 - sqr(3.79E-4)*beta_*c*mag(Enorm)*1E-2/sqr(wf_) ; // 1E-2 is for converting V/m to V/cm

        //Info << "vofy = " << vofy << endl;

        Gamma_FE = c*1.54E-6/1.602e-19*sqr(beta_*c*mag(Enorm) )/1.1/wf_*exp(-6.85E9*pow(wf_,1.5)*vofy/beta_/(c*mag(Enorm) + SMALL) ) ;

        //Info << "Gamma_FE = " << Gamma_FE << endl;
    }

    //vectorField temp = 0.0*n;

    vectorField temp = (0.25*sqrt(8.0*1.38e-23*Tef/9.1e-31/acos(-1.0)) - (Gamma_se + Gamma_FE)/Nef)*n;

    operator == (temp);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::electronThermalVelocity::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
    os.writeKeyword("Edepend")
        << Edepend_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_emission")
        << FE_ << token::END_STATEMENT << nl;
    os.writeKeyword("field_enhancement_factor")
        << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("work_function")
        << wf_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        electronThermalVelocity
    );
}


// ************************************************************************* //
