# SOMAFOAM

This software consists of FOAM based FV method and other utilities designed and implemented for modular, multiphysics plasma fluid simulation.

The foam-extend snapshot used in this project corresponds to following,
```
commit efcc2b1b7df8543c7873f89a6f50c30e047b1b11
Date:   Thu Apr 12 14:05:28 2018 +0100
```
There were certain changes made to the foam-extend base depending on our needs.
The code compiles on newer gcc versions (including 10.2.1)

## Building

The repository contains sources only; compiled libraries and applications are
machine-specific (compiler, glibc) and are built locally.

Requirements: gcc/g++, make, flex and OpenMPI (with development headers).
Tested on Debian 11 (gcc 10.2.1, OpenMPI 4.1.0).

```
git clone https://github.com/vayyaswamy/SOMAFOAM
cd SOMAFOAM
./install.sh            # builds wmake, all libraries and applications
```

The libraries go to `lib/linux64GccDPOpt` and the applications to
`bin/linux64GccDPOpt`. Before running, load the environment in each new shell:

```
source <path-to>/SOMAFOAM/etc/bashrc
```

Run the cases in `examples/plasma` with `somaFoam` (some `controlDict` files
still name the older `plasmaSimFoam`, which no longer runs them).

## Inner iterations and field relaxation

The number of passes per time step over each species equation and the
electron temperature equation, which resolve the nonlinear chemistry source,
can be set in `constant/plasmaProperties` (defaults shown). A loop that ends
above the tolerance is reported in the log.

```
innerIterations
{
    chargedSpecies      4;
    neutralSpecies      2;
    electronTemperature 6;
    tolerance           1e-5;   // initial residual that ends the passes
    reportUnconverged   yes;
}
```

Field relaxation factors in `system/fvSolution` follow the usual convention:
`<field>Final` applies on the final PIMPLE iteration of a time step, so that

```
fields { ".*" 0.8; ".*Final" 1.0; }
```

relaxes intermediate iterations only and a run with a single PIMPLE iteration
is not relaxed. (Before this was honoured, every time step took 80 % of its
change with these settings.)

## Acceleration of slow neutral species

Metastables reach their periodic steady state over times that are thousands
of periods of the applied voltage. They can be advanced alone, with a large
time step, between blocks of full simulation (`constant/plasmaProperties`):

```
slowSpeciesAcceleration
{
    species         (Arm);
    period          2.5e-8;   // period of the applied voltage [s]
    fullCycles      20;       // periods of full simulation per block
    averageCycles   2;        // periods at the end of a block over which
                              // the chemistry source is averaged
    deltaT          1e-6;     // time step of the advance [s]
    nSteps          20;       // steps per advance
    maxChangeFactor 2;        // limit on the density change per advance
    tolerance       1e-3;     // advances stop below this relative change
}
```

After each block the listed species are advanced by `nSteps*deltaT` with the
period-averaged source, linearised in their own density, and diffusion; the
charged species, the electron temperature and the potential are left
unchanged and adjust during the next block. Only neutral species transported
by diffusion can be listed. The time reported by the solver does not include
the advances.

## Electrode voltage and current

The `electrodeVoltageCurrent` function object writes, for each listed patch,
the voltage and the conduction, displacement and total current versus time to
`<case>/<name>/<start time>/<patch>.dat` (currents positive from the electrode
into the plasma; the last column is the total current density). In
`system/controlDict`:

```
functions
{
    electrodes
    {
        type                electrodeVoltageCurrent;
        functionObjectLibs  ("libplasmaFunctionObjects.so");
        patches             (electrode ground);
        // optional: outputInterval 1;
    }
}
```

## Adaptive mesh refinement (1D)

`somaFoam` can refine and coarsen a one-dimensional mesh during the run. It is
switched on per case by adding `constant/dynamicMeshDict`; without that file
the mesh is static and results are unchanged.

```
dynamicFvMesh   dynamicRefine1DFvMesh;

dynamicRefine1DFvMeshCoeffs
{
    direction           (1 0 0);   // direction of the 1D mesh
    refineInterval      50;        // time steps between mesh updates
    indicators                     // indicator = maximum over the entries
    (
        { type relativeGradient; field N_electron; floor 1e12; }
        { type relativeGradient; field N_Arp1;     floor 1e12; }
    );
    weightedAverages    ((Te N_electron));
    lowerRefineLevel    0.2;       // refine where indicator is above this
    upperRefineLevel    1e30;
    unrefineLevel       0.07;      // merge where indicator is below this
    nBufferLayers       2;
    maxRefinement       2;         // each level halves the cell size
    maxCells            2000;
}
```

and loading the library in `system/controlDict`:

```
libs ( "libfoam.so" "liblduSolvers.so" "libplasmaCookBook.so" "libplasmaDynamicMesh.so" );
```

Indicator types, for any `volScalarField` or `volVectorField` of the solver
(`N_<specie>`, `Te`, `Phi`, `E`, ...), each with an optional `weight`:

| type | value |
|---|---|
| `relativeGradient` | `|f1 - f2| / (0.5(|f1| + |f2|) + floor)` between neighbouring cells, i.e. the relative gradient times the cell size |
| `gradient` | `|f1 - f2| / scale` |
| `magnitude` | `|f| / scale` |

Keep `floor` well below the densities in the sheaths, otherwise they are not
refined, and `unrefineLevel` below half of `lowerRefineLevel`.

Limits: 1D meshes only (a single row of hexahedral cells; other meshes stop
with an error), `solutionDomain plasma` only (not with dielectric regions),
and the `temporal` chemistry mode has not been adapted. Tested with the
`driftDiffusion` model on the 100 Torr argon case at 100 V: 75 cells with two
levels reproduce a uniform 300-cell mesh within 3 %.

## Adaptive mesh refinement (2D and 3D)

For 2D (one cell thick, `empty` front and back) and 3D meshes `somaFoam` uses
the polyhedral refinement of foam-extend 4.1 (`dynamicPolyRefinementFvMesh`,
ported into `src/dynamicMesh`): cells are split in four (2D) or eight (3D)
with hanging nodes, and merged again. The same indicators as in 1D are
available through the `plasmaIndicatorRefinement` selection.
`constant/dynamicMeshDict`:

```
dynamicFvMesh   dynamicPolyRefinementFvMesh;

dynamicPolyRefinementFvMeshCoeffs
{
    refineInterval      50;        // time steps between refinements
    unrefineInterval    50;
    separateUpdates     false;

    active              yes;
    maxCells            200000;
    maxRefinementLevel  2;
    nRefinementBufferLayers   2;
    nUnrefinementBufferLayers 4;
    edgeBasedConsistency yes;

    refinementSelection
    {
        type            plasmaIndicatorRefinement;
        indicators
        (
            { type relativeGradient; field N_electron; floor 1e12; }
            { type relativeGradient; field N_Arp1;     floor 1e12; }
        );
        lowerRefineLevel    0.2;
        unrefineLevel       0.07;
    }
}
```

and in `system/controlDict`:

```
libs ( "libfoam.so" "liblduSolvers.so" "libplasmaCookBook.so" "libtopoChangerFvMesh.so" "libplasmaDynamicMesh.so" );
```

The 4.1 selections (`fieldBoundsRefinement`, `minCellSizeRefinement`,
`compositeRefinementSelection`, ...) can be used as well.

With hanging nodes the faces between fine and coarse cells are non-orthogonal:
use `corrected` Laplacian and `snGrad` schemes in `system/fvSchemes`.

`examples/plasma/2DAdaptiveMeshArgon` is a worked example (seeded plasma blob,
one refinement level).

Parallel runs keep their decomposition while they run (the 4.1 load balancing
is not ported), so refinement can leave the processors unevenly loaded. To
balance them, stop the run, call `rebalancePar` in the case directory and
restart:

```
mpirun -np 4 somaFoam -parallel     # stops at endTime
rebalancePar                        # needs startFrom latestTime in controlDict
mpirun -np 4 somaFoam -parallel     # after raising endTime
```

`rebalancePar` rebuilds the refined mesh and the fields from the processor
directories (`reconstructParMesh`), decomposes them again according to
`system/decomposeParDict` (the number of processors may be changed) and
transfers the refinement levels with the `refinementLevelsPar` utility. The
old processor directories are kept in `beforeRebalance_<time>`. Refined cells
whose siblings end up on different processors cannot be merged again, so the
mesh may stay slightly finer along the new processor boundaries.

Limits: `solutionDomain plasma` only; merged cells take the volume average of
all fields.
Tested with the `driftDiffusion` model: the 100 V argon case on a thin 2D
mesh matches the 1D result within 1 %, the seeded-blob example matches a
uniform fine mesh within 2 % (0.4 % RMS), and runs restart from a refined
mesh. On 2 and 4 processors the refined mesh is identical to the
serial one and the fields agree with the serial run to 5e-5 of their peak
(5e-4 for the electron temperature).

Contributors:
1) Venkattraman Ayyaswamy (https://me.ucmerced.edu/content/venkattraman-venkatt-ayyaswamy)
2) Abhishek Kumar Verma 
3) Saurav Gautam (https://me.ucmerced.edu/content/saurav-gautam)
4) Jose Alfredo Millan Higuera (https://www.linkedin.com/in/jmillanhiguera/)

For research conducted using SOMAFOAM please cite this paper: https://doi.org/10.1016/j.cpc.2021.107855

Paper that compares SOMAFOAM results with experiment: 
https://doi.org/10.1063/5.0041386
