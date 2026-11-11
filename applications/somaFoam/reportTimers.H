// Report of the wall-clock times of solverTimers.H
{
    const double total = solverTimers::now() - solverTimers::start;
    const double chem = Foam::multiSpeciesPlasmaModel::chemistryWallTime;
    const double other =
        total - solverTimers::poisson - solverTimers::species
      - solverTimers::energy - solverTimers::surface;

    Info<< "Timing (wall s): total " << total
        << ", Poisson " << solverTimers::poisson
        << ", plasma species " << solverTimers::species
        << ", electron energy " << solverTimers::energy
        << ", surface charge " << solverTimers::surface
        << ", other " << other
        << "; chemistry (in species and energy) " << chem
        << "; outer iterations " << solverTimers::outer
        << " in " << solverTimers::steps << " steps" << endl;
}
