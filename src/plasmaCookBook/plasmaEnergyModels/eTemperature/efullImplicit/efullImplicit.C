/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "efullImplicit.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(efullImplicit, 0);
    addToRunTimeSelectionTable(eTemp, efullImplicit, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efullImplicit::efullImplicit
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

Foam::scalar Foam::efullImplicit::correct
(
    psiChemistryModel& chemistry,
    const volVectorField& E
)
{

    lduSolverPerformance solverPerf;
    volScalarField& TeC = thermo().Te();
    TeC.storePrevIter();
    scalar initialResidual = 1.0;

    int icorr = 0;

    while
    (
        (icorr++ < mspm().nCorrTe())
     && (initialResidual >= mspm().innerTolerance())
    )

    {

       eeFlux = 2.5*plasmaConstants::boltzC*mspm().J(eIndex_);

       surfaceScalarField eeFluxF = fvc::interpolate(eeFlux) & mesh().Sf();

       // Electron flux times electric field
       volScalarField jDotE = mspm().J(eIndex_) & E;

       if (mspm().hasFaceFlux(eIndex_))
       {
           // The flux of the electron equation through the faces carries
           // the energy and does the work against the field:
           // J.E = -J.grad(Phi) = Phi div(J) - div(J Phi)
           const surfaceScalarField& Jf = mspm().faceFlux(eIndex_);

           const volScalarField& Phi =
               mesh().lookupObject<volScalarField>("Phi");

           eeFluxF = 2.5*plasmaConstants::boltzC*Jf;

           // Potential at the face centres. Where the line between two
           // cell centres does not pass through the face centre (faces
           // between fine and coarse cells of a refined mesh) the
           // interpolated value is corrected with the gradient, -E
           surfaceScalarField PhiF(linear<scalar>(mesh()).interpolate(Phi));

           const skewCorrectionVectors& skew =
               skewCorrectionVectors::New(mesh());

           if (skew.skew())
           {
               PhiF -= skew() & linear<vector>(mesh()).interpolate(E);
           }

           jDotE.internalField() =
               Phi.internalField()*fvc::div(Jf)().internalField()
             - fvc::div(Jf*PhiF)().internalField();
       }

       volScalarField eeSource = - plasmaConstants::eCharge*jDotE - mspm().electronTempSource(chemistry);

        volScalarField eeSource_Su = plasmaConstants::eCharge*jDotE
                                    + mspm().electronTempSource(chemistry)
                                    - mspm().dElectronTempSourceDTe(chemistry)*TeC;

        volScalarField eeSource_SuSp = mspm().dElectronTempSourceDTe(chemistry);


       const volScalarField& Ne = mspm().N(eIndex_);

        const volScalarField eConductivity
        (
            mspm().electronConductivity(chemistry)
        );

        // Convection with the electron flux and conduction
        tmp<fvScalarMatrix> transport
        (
            mspm().scharfetterGummel(eIndex_)
          ? scharfetterGummel(eeFluxF, eConductivity, TeC)
          : fvm::div(eeFluxF, TeC)
          - fvm::laplacian(eConductivity, TeC, "laplacian(eC,Te)")
        );

        fvScalarMatrix TeEqn
        (
            fvm::ddt((1.5*plasmaConstants::boltzC*Ne), TeC)
            + transport
            + fvm::SuSp(eeSource_SuSp, TeC)
            + eeSource_Su
        );

       solverPerf = TeEqn.solve();

        initialResidual = solverPerf.initialResidual();

        mspm().limitField(TeC, mspm().TeMin(), mspm().TeMax());

        thermo().correct();

        mspm().updateTemperature();

        mspm().updateChemistryCollFreq(chemistry);

    }

    mspm().reportInnerIterations("Te", mspm().nCorrTe(), initialResidual);

    // Solved once per time step, after the PIMPLE loop
    relaxFinal(TeC, true);

    mspm().limitField(TeC, mspm().TeMin(), mspm().TeMax());

    thermo().correct();

    mspm().updateTemperature();

    mspm().updateChemistryCollFreq(chemistry);

    volScalarField divEeFlux = mag(0.5*fvc::div(eeFlux));

    scalar maxDivEeFlux = gMax(divEeFlux);

    return maxDivEeFlux;
}


// ************************************************************************* //
