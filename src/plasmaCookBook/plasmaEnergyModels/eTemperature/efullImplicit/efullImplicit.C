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

    // The species update leaves the rates of the last species unchanged:
    // the energy source of the reactions for the present state
    if (!mspm().chemistryPerPass())
    {
        mspm().updateElectronEnergyChemistry(chemistry);
    }

    // Electron energy density Ne*Te before each pass: a pass ends the
    // iterations when it changes Ne*Te by less than innerTolerance relative to
    // its largest value (cells without electrons, where Te is not determined,
    // do not hold it up; the normalised residual used before did)
    scalarField NeTePass
    (
        mspm().N(eIndex_).internalField()*thermo().Te().internalField()
    );

    const dimensionedScalar rDeltaT(1.0/mesh().time().deltaT());

    // Energy density of the start of the step over the step
    const volScalarField oldEnergy
    (
        rDeltaT*1.5*plasmaConstants::boltzC
       *mspm().N(eIndex_).oldTime()*TeC.oldTime()
    );

    // The energy flux of the electrons and the work of the field depend
    // only on the electron flux and the potential, which the passes do not
    // change
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

           // The face flux of the electrons has no part carried by the gas
           // flow (with the cell flux J, addGasFlux adds it)
           if (mesh().foundObject<surfaceScalarField>("phiGas"))
           {
               eeFluxF +=
                   2.5*plasmaConstants::boltzC
                  *mesh().lookupObject<surfaceScalarField>("phiGas")
                  *fvc::interpolate(mspm().N(eIndex_));
           }

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


    while
    (
        (icorr++ < mspm().nCorrTe())
     && (initialResidual >= mspm().innerTolerance())
    )

    {

        // Energy source of the reactions and its Te derivative, once each
        const volScalarField eeChemSource
        (
            mspm().electronTempSource(chemistry)
        );

        volScalarField eeSource_SuSp = mspm().dElectronTempSourceDTe(chemistry);

        volScalarField eeSource_Su = plasmaConstants::eCharge*jDotE
                                    + eeChemSource
                                    - eeSource_SuSp*TeC;


       const volScalarField& Ne = mspm().N(eIndex_);

        const volScalarField eConductivity
        (
            mspm().electronConductivity(chemistry)
        );

        // The sources keep the Te of the start of the pass through the
        // non-orthogonal correctors
        for (label nonOrth = 0; nonOrth <= mspm().nNonOrthCorr(); nonOrth++)
        {
            // Convection with the electron flux and conduction
            tmp<fvScalarMatrix> transport
            (
                mspm().scharfetterGummel(eIndex_)
              ? scharfetterGummel(eeFluxF, eConductivity, TeC)
              : fvm::div(eeFluxF, TeC)
              - fvm::laplacian(eConductivity, TeC, "laplacian(eC,Te)")
            );

            // Time derivative of the electron energy density 3/2 k Ne Te,
            // (Ne Te - Ne^0 Te^0)/dt. fvm::ddt(rho, Te) would take the old
            // time of rho from the temporary 3/2 k Ne, a copy of the new
            // value, and so discretize 3/2 k Ne dTe/dt, dropping
            // 3/2 k Te dNe/dt
            fvScalarMatrix TeEqn
            (
                fvm::Sp(rDeltaT*1.5*plasmaConstants::boltzC*Ne, TeC)
              - oldEnergy
              + transport
              + fvm::SuSp(eeSource_SuSp, TeC)
              + eeSource_Su
            );

            solverPerf = TeEqn.solve();

            mspm().limitField(TeC, mspm().TeMin(), mspm().TeMax());
        }

        initialResidual = mspm().passChange
        (
            mspm().N(eIndex_).internalField()*TeC.internalField(),
            NeTePass
        );

        // The thermo (gas temperature, compressibility, viscosity) does not
        // depend on Te and is not updated here

        mspm().updateTemperature();

        // The passes need only the electron energy source; all rates are
        // updated after the loop
        mspm().updateElectronEnergyChemistry(chemistry);
    }

    mspm().reportInnerIterations("Te", mspm().nCorrTe(), initialResidual);

    // Solved once per time step, after the PIMPLE loop
    relaxFinal(TeC, true);

    mspm().limitField(TeC, mspm().TeMin(), mspm().TeMax());

    mspm().updateTemperature();

    // Terms of the energy equation at the solution, as gains [W/m3], with
    // "writeEnergyBudget yes;" in controlDict at the write times. The
    // residual is the imbalance: what the limits on Te and the
    // linearisation of the reaction source left
    if
    (
        mesh().time().outputTime()
     && mesh().time().controlDict().lookupOrDefault<Switch>
        (
            "writeEnergyBudget", false
        )
    )
    {
        const volScalarField& Ne = mspm().N(eIndex_);

        const volScalarField eConductivity
        (
            mspm().electronConductivity(chemistry)
        );

        volScalarField eStorage
        (
            "eStorage",
            rDeltaT*1.5*plasmaConstants::boltzC*Ne*TeC - oldEnergy
        );
        volScalarField eHeating("eHeating", -plasmaConstants::eCharge*jDotE);
        volScalarField eChem("eChem", -mspm().electronTempSource(chemistry));
        volScalarField eConvection("eConvection", -fvc::div(eeFluxF, TeC));
        volScalarField eConduction
        (
            "eConduction",
            fvc::laplacian(eConductivity, TeC, "laplacian(eC,Te)")
        );
        volScalarField eResidual
        (
            "eResidual",
            eStorage - eHeating - eChem - eConvection - eConduction
        );

        eStorage.write();
        eHeating.write();
        eChem.write();
        eConvection.write();
        eConduction.write();
        eResidual.write();
    }

    mspm().updateChemistryCollFreq(chemistry);

    mspm().setChemistryCurrent();

    volScalarField divEeFlux = mag(0.5*fvc::div(eeFlux));

    scalar maxDivEeFlux = gMax(divEeFlux);

    return maxDivEeFlux;
}


// ************************************************************************* //
