// Step controls of the solutionDomain plasma branch with a converged outer
// loop (plasmaConvergedStep.H); the same as those of the plasmaDielectric
// branch (somaFoam.C), which declares its own.

// Upper limit of the step after a rejected step; GREAT when none
scalar deltaTRejectionCap = GREAT;

// Optional bound on the ionization number nu dt (controlDict
// maxIonizationNumber, default 0 = off), over the cells with at least
// ionizationNumberMinDensity electrons and Te below TeMax
const scalar maxIonizationNumber =
    runTime.controlDict().lookupOrDefault<scalar>("maxIonizationNumber", 0);

const scalar ionizationNumberMinDensity =
    runTime.controlDict().lookupOrDefault<scalar>
    (
        "ionizationNumberMinDensity",
        0
    );

const label ionizationElectronI = mspm().species()["electron"];

// Largest positive net ionization frequency of the electrons [1/s]
scalar nuIonizationMax = 0;

// Optional upper limit of the step while the discharge conducts
// (controlDict deltaTMaxDischarge [s], dischargeCurrent [A])
const scalar deltaTMaxDischarge =
    runTime.controlDict().lookupOrDefault<scalar>("deltaTMaxDischarge", GREAT);

const scalar dischargeCurrent =
    runTime.controlDict().lookupOrDefault<scalar>("dischargeCurrent", 0);
