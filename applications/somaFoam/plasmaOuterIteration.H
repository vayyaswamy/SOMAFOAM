// One outer iteration of the solutionDomain plasma step: Poisson equation,
// plasma equations and the surface charge of the dielectric side walls

solverTimers::outer++;

{
    double tTimer = solverTimers::now();

    while (pimple.correct())
    {
        #include "solvePoisson.H"
    }

    solverTimers::poisson += solverTimers::now() - tTimer;
    tTimer = solverTimers::now();

    #include "plasmaEqn.H"

    solverTimers::species += solverTimers::now() - tTimer;
    tTimer = solverTimers::now();

    #include "surfaceCharge_new.H"

    solverTimers::surface += solverTimers::now() - tTimer;
}
