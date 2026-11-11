// Wall-clock time of the parts of a step, reported with "timingReport yes;"
// in system/controlDict (at each write and at the end of the run). The
// chemistry time is accumulated by the plasma model
// (multiSpeciesPlasmaModel::chemistryWallTime) and is part of the species
// and electron energy times.

#ifndef solverTimers_H
#define solverTimers_H

#include <chrono>

namespace solverTimers
{
    inline double now()
    {
        return std::chrono::duration<double>
        (
            std::chrono::steady_clock::now().time_since_epoch()
        ).count();
    }

    static double poisson = 0;
    static double species = 0;
    static double energy = 0;
    static double surface = 0;
    static double start = 0;
    static long outer = 0;
    static long steps = 0;
}

#endif
