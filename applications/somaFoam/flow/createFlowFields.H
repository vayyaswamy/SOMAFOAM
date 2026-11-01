// Gas flow (optional), solved with the equations of the compressible flow
// solvers of foam-extend 4.0: rhoSimpleFoam for a steady flow and
// rhoPimpleFoam for a transient one, with the turbulence models of
// src/turbulenceModels/compressible (constant/turbulenceProperties: laminar,
// RAS or LES). The plasma enters only through the force e sum_i z_i N_i E on
// the space charge and the heating of the gas (ion Joule heating and
// ion-gas energy exchange, multiSpeciesPlasmaModel::ionTempSource), and the
// species are carried by the gas (phiGas, multiSpeciesPlasmaModel). The gas
// energy is the sensible enthalpy hs of the plasma thermo.
//
// constant/flowProperties (without the file, or with flow off, nothing here
// is solved):
//
//     flow            on;
//     solution        steady; // steady (default), transient or frozen
//     energy          on;     // gas energy, heated by the plasma
//     bodyForce       (0 0 0); // uniform force [N/m3] (optional)
//     heatSource      0;      // uniform heating [W/m3] (optional)
//
//     // steady: solved at the start by SIMPLE iterations (SIMPLE in
//     // fvSolution, with rhoMin, rhoMax and the relaxation factors of
//     // rhoSimpleFoam)
//     steadyMaxIterations 2000;
//     steadyTolerance 1e-6;   // largest relative change of U, p (and T)
//                             // in an iteration
//     updateInterval  0;      // plasma time [s] after which the steady flow
//                             // is solved again with the plasma force and
//                             // heating averaged over it; 0: never
//
//     // transient: one rhoPimpleFoam step per plasma step (flowPIMPLE in
//     // fvSolution, with pMin and the controls of PIMPLE)
//
// A steady flow at the start (without plasma force and heating) is the flow
// set by the boundaries, such as a gas feed through a reactor; the plasma
// then sees a fixed flow. With updateInterval it follows the period-averaged
// plasma (electrohydrodynamic flow, gas heating) on the slow time scale of
// the gas. frozen holds the flow of 0/U and 0/p.
//
// The fields and dictionaries are those of rhoSimpleFoam and rhoPimpleFoam:
// 0/U, constant/turbulenceProperties (and RASProperties or LESProperties),
// the solvers and schemes of U, p and hs. With energy on, the gas temperature
// model of plasmaProperties should be none.

IOdictionary flowProperties
(
    IOobject
    (
        "flowProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const Switch flowOn(flowProperties.lookupOrDefault<Switch>("flow", false));

const word flowSolution
(
    flowProperties.lookupOrDefault<word>("solution", "steady")
);

if
(
    flowOn
 && flowSolution != "steady"
 && flowSolution != "transient"
 && flowSolution != "frozen"
)
{
    FatalErrorIn(args.executable())
        << "flowProperties: unknown solution " << flowSolution
        << " (steady, transient or frozen)" << exit(FatalError);
}

// Without a gas flow and with the gas energy model none, nothing changes the
// gas: its thermo (temperature, density) need not be updated after the
// species (plasmaEqn.H)
const bool gasThermoFixed =
    !flowOn
 && plasmaProperties.subDict("energyModel").lookupOrDefault<word>("gas", "none")
 == "none";

const bool flowSteady = flowOn && flowSolution == "steady";
const bool flowTransient = flowOn && flowSolution == "transient";

const Switch flowEnergy
(
    flowProperties.lookupOrDefault<Switch>("energy", true)
);

const label steadyMaxIterations
(
    flowProperties.lookupOrDefault<label>("steadyMaxIterations", 2000)
);

const scalar steadyTolerance
(
    flowProperties.lookupOrDefault<scalar>("steadyTolerance", 1e-6)
);

const scalar flowUpdateInterval
(
    flowProperties.lookupOrDefault<scalar>("updateInterval", 0)
);

const dimensionedVector flowBodyForce
(
    "bodyForce",
    dimForce/dimVolume,
    flowProperties.lookupOrDefault<vector>("bodyForce", vector::zero)
);

const dimensionedScalar flowHeatSource
(
    "heatSource",
    dimEnergy/dimTime/dimVolume,
    flowProperties.lookupOrDefault<scalar>("heatSource", 0)
);

autoPtr<volVectorField> UgasPtr;
autoPtr<surfaceScalarField> phiPtr;
autoPtr<surfaceScalarField> phiGasPtr;
autoPtr<compressible::turbulenceModel> turbulence;
autoPtr<volScalarField> DpDtPtr;
autoPtr<simpleControl> flowSimple;
autoPtr<pimpleControl> flowPimple;

// Controls of rhoSimpleFoam (steady) and rhoPimpleFoam (transient)
label pRefCell = 0;
scalar pRefValue = 0.0;
dimensionedScalar rhoMin("rhoMin", dimDensity, 0);
dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT);
dimensionedScalar pMin("pMin", dimPressure, 0);
dimensionedScalar initialMass("initialMass", dimMass, 0);

#include "initContinuityErrs.H"

// Plasma force and heating, averaged over the plasma steps since the last
// steady solution (updateInterval), or of the last plasma step (transient)
autoPtr<volVectorField> plasmaForcePtr;
autoPtr<volScalarField> plasmaHeatPtr;
autoPtr<vectorField> plasmaForceSum;
autoPtr<scalarField> plasmaHeatSum;
scalar flowTimeSum = 0;

if (flowOn)
{
    Info<< "Gas flow: " << flowSolution;

    if (flowSteady && flowUpdateInterval > 0)
    {
        Info<< ", solved again every " << flowUpdateInterval << " s";
    }

    Info<< ", gas energy " << flowEnergy << endl;

    UgasPtr.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    phiPtr.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rho*UgasPtr()) & mesh.Sf()
        )
    );

    // Volume flux of the gas, which carries the species (looked up by the
    // plasma model)
    phiGasPtr.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiGas",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiPtr()/linearInterpolate(rho)
        )
    );

    mesh.schemesDict().setFluxRequired(thermo.p().name());

    plasmaForcePtr.reset
    (
        new volVectorField
        (
            IOobject
            (
                "plasmaForce",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimVolume, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    plasmaHeatPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "plasmaHeat",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    plasmaForceSum.reset(new vectorField(mesh.nCells(), vector::zero));
    plasmaHeatSum.reset(new scalarField(mesh.nCells(), 0.0));

    if (flowSolution != "frozen")
    {
        Info<< "Creating turbulence model\n" << endl;

        turbulence = compressible::turbulenceModel::New
        (
            rho,
            UgasPtr(),
            phiPtr(),
            thermo
        );
    }

    if (flowSteady)
    {
        flowSimple.reset(new simpleControl(mesh));

        setRefCell(thermo.p(), flowSimple().dict(), pRefCell, pRefValue);

        rhoMax = dimensionedScalar(flowSimple().dict().lookup("rhoMax"));
        rhoMin = dimensionedScalar(flowSimple().dict().lookup("rhoMin"));

        initialMass = fvc::domainIntegrate(rho);
    }

    if (flowTransient)
    {
        flowPimple.reset(new pimpleControl(mesh, "flowPIMPLE"));

        pMin = dimensionedScalar(flowPimple().dict().lookup("pMin"));

        DpDtPtr.reset
        (
            new volScalarField
            (
                "DpDt",
                fvc::DDt
                (
                    surfaceScalarField("phiU", phiPtr()/fvc::interpolate(rho)),
                    thermo.p()
                )
            )
        );
    }

    if (flowEnergy)
    {
        // The heat capacity of thermoProperties is per kg; a molar value per
        // mol (about 20.8 for a monatomic gas) gives a gas that heats about
        // 25 times too fast
        const scalar CpMin = gMin(thermo.Cp()().internalField());

        if (CpMin < 100)
        {
            WarningIn(args.executable())
                << "Gas heat capacity " << CpMin << " J/(kg K): the heat"
                << " capacities in thermoProperties are per kg"
                << " (argon: 520 J/(kg K))." << endl;
        }

        const word gasModel
        (
            plasmaProperties.subDict("energyModel")
           .lookupOrDefault<word>("gas", "none")
        );

        if (gasModel != "none")
        {
            WarningIn(args.executable())
                << "Gas energy is solved by the flow (flowProperties energy"
                << " on) and by the gas temperature model " << gasModel
                << " of plasmaProperties." << nl
                << "    Set gas none in energyModel." << endl;
        }
    }

    // Steady flow at the start, without the plasma
    if (flowSteady)
    {
        #include "flowSteady.H"
    }
}
