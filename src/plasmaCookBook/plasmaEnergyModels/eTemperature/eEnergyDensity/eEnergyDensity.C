/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "eEnergyDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
#include "skewCorrectionVectors.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eEnergyDensity, 0);
    addToRunTimeSelectionTable(eTemp, eEnergyDensity, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eEnergyDensity::eEnergyDensity
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E,
    const dictionary& dict
)
:
    eTemp(thermo, mspm, E),
    eeFlux
    (
        IOobject
        (
            "eeFlux",
            thermo.T().mesh().time().timeName(),
            thermo.T().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.T().mesh(),
        dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
    ),
    eSpecie("electron"),
    eIndex_(mspm.species()[eSpecie])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::eEnergyDensity::correct
(
    psiChemistryModel& chemistry,
    const volVectorField& E
)
{
    volScalarField& TeC = thermo().Te();
    const volScalarField& Ne = mspm().N(eIndex_);
    const dimensionedScalar k15(1.5*plasmaConstants::boltzC);

    // The species update leaves the rates of the last species unchanged:
    // the energy source of the reactions for the present state
    if (!mspm().chemistryPerPass())
    {
        mspm().updateElectronEnergyChemistry(chemistry);
    }

    // Energy density before each pass (see passChange)
    scalarField NeTePass(Ne.internalField()*TeC.internalField());

    const dimensionedScalar rDeltaT(1.0/mesh().time().deltaT());

    // Energy density of the start of the step over the step
    const volScalarField oldEnergy
    (
        rDeltaT*k15*Ne.oldTime()*TeC.oldTime()
    );

    // Drift velocity and diffusivity of the energy, 5/3 of those of the
    // electrons
    const volVectorField& Fe =
        mesh().lookupObject<volVectorField>("F_" + eSpecie);
    const volScalarField Deps
    (
        (5.0/3.0)*mesh().lookupObject<volScalarField>("D_" + eSpecie)
    );

    // Work of the field on the electrons, with the face flux of the electron
    // equation where there is one: J.E = Phi div(J) - div(J Phi)
    volScalarField jDotE(mspm().J(eIndex_) & E);

    if (mspm().hasFaceFlux(eIndex_))
    {
        const surfaceScalarField& Jf = mspm().faceFlux(eIndex_);

        const volScalarField& Phi = mesh().lookupObject<volScalarField>("Phi");

        surfaceScalarField PhiF(linear<scalar>(mesh()).interpolate(Phi));

        const skewCorrectionVectors& skew = skewCorrectionVectors::New(mesh());

        if (skew.skew())
        {
            PhiF -= skew() & linear<vector>(mesh()).interpolate(E);
        }

        jDotE.internalField() =
            Phi.internalField()*fvc::div(Jf)().internalField()
          - fvc::div(Jf*PhiF)().internalField();
    }

    // Energy density: the wall value from the Te boundary condition, coupled
    // and empty patches as for Te
    wordList patchTypes(TeC.boundaryField().size());

    forAll(patchTypes, patchi)
    {
        const fvPatchScalarField& pTe = TeC.boundaryField()[patchi];

        patchTypes[patchi] =
            pTe.coupled() || pTe.type() == "empty"
          ? pTe.type()
          : word(fixedValueFvPatchScalarField::typeName);
    }

    volScalarField nEps
    (
        IOobject
        (
            "electronEnergyDensity",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", k15.dimensions()*Ne.dimensions()*TeC.dimensions(), 0),
        patchTypes
    );

    nEps.internalField() = k15.value()*Ne.internalField()*TeC.internalField();

    scalar initialResidual = 1.0;
    int icorr = 0;

    while
    (
        (icorr++ < mspm().nCorrTe())
     && (initialResidual >= mspm().innerTolerance())
    )
    {
        // Wall values: 3/2 kB n_e Te with the wall Te of its condition
        TeC.correctBoundaryConditions();

        forAll(nEps.boundaryField(), patchi)
        {
            fvPatchScalarField& pEps = nEps.boundaryField()[patchi];

            if (!pEps.coupled() && pEps.type() != "empty")
            {
                pEps ==
                    k15.value()*Ne.boundaryField()[patchi]
                   *TeC.boundaryField()[patchi];
            }
        }

        nEps.correctBoundaryConditions();

        // Energy loss of the reactions and its Te derivative; dS/dn_eps =
        // dS/dTe/(3/2 kB n_e)
        const volScalarField eeChemSource(mspm().electronTempSource(chemistry));

        const volScalarField dSdTe(mspm().dElectronTempSourceDTe(chemistry));

        const volScalarField rk15Ne
        (
            1.0/(k15*max(Ne, dimensionedScalar("small", Ne.dimensions(), VSMALL)))
        );

        // Te of the source linearisation is that of the start of the pass,
        // through the non-orthogonal correctors
        for (label nonOrth = 0; nonOrth <= mspm().nNonOrthCorr(); nonOrth++)
        {
            fvScalarMatrix epsEqn
            (
                fvm::Sp(rDeltaT, nEps) - oldEnergy
              + mspm().driftDiffusionTerms(eIndex_, Fe, Deps, nEps, 5.0/3.0)
              + mspm().gasConvection(nEps, "div(phiGas,Ni)")
              + fvm::SuSp(dSdTe*rk15Ne, nEps)
              + plasmaConstants::eCharge*jDotE
              + eeChemSource
              - dSdTe*TeC
            );

            epsEqn.solve(mesh().solutionDict().solver(TeC.name()));
        }

        // Te from the energy density, within its limits
        TeC.internalField() = nEps.internalField()*rk15Ne.internalField();

        mspm().limitField(TeC, mspm().TeMin(), mspm().TeMax());

        nEps.internalField() =
            k15.value()*Ne.internalField()*TeC.internalField();

        initialResidual = mspm().passChange
        (
            Ne.internalField()*TeC.internalField(),
            NeTePass
        );

        mspm().updateTemperature();

        // The passes need only the electron energy source; all rates are
        // updated after the loop
        mspm().updateElectronEnergyChemistry(chemistry);
    }

    mspm().reportInnerIterations("Te", mspm().nCorrTe(), initialResidual);

    TeC.correctBoundaryConditions();

    mspm().updateTemperature();

    mspm().updateChemistryCollFreq(chemistry);

    mspm().setChemistryCurrent();

    // The time step limit of efullImplicit, from 5/2 kB J_e
    eeFlux = 2.5*plasmaConstants::boltzC*mspm().J(eIndex_);

    return gMax(mag(0.5*fvc::div(eeFlux))().internalField());
}


// ************************************************************************* //
