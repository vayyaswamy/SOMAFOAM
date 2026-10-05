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

Run the cases in `examples/plasma` with `somaFoam`.

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

## Plasma with dielectric regions

Two solvers handle a plasma bounded by dielectrics (`solutionDomain
plasmaDielectric` in `constant/electroMagnetics`, regions listed in
`constant/regionProperties`, interfaces as `regionCouple` patches). They use
the same case and the same plasma step and differ in how the potential is
coupled across the interfaces:

| Solver | Coupling | Potential on the interface patches |
|---|---|---|
| `somaFoam` | plasma and dielectrics in one matrix | `coupledPotential` |
| `plasmaMultiRegionFoam` | regions solved in turn and iterated within each time step | `iterativeCoupledPotential` |

For `plasmaMultiRegionFoam`, set the type of `Phi` on both sides of every
interface (`0/Phi` and `0/<dielectric>/Phi`) to `iterativeCoupledPotential`,
keeping the `remoteField` and `surfaceCharge` entries, and provide a linear
solver for `Phi` in `system/fvSolution` of every region. The plasma side
fixes the interface potential and each dielectric imposes the normal gradient
that satisfies Gauss's law with the surface charge; the interface potential
is relaxed (Aitken) until both agree. Optional controls in the plasma's
`system/fvSolution`:

```
plasmaDielectricCoupling
{
    maxIterations       50;
    tolerance           1e-8;   // change of the interface potential relative
                                // to the largest potential
    initialRelaxation   0.5;
}
```

On `examples/plasmaDielectric/ArgonDBD` the iteration takes 3 to 5 passes per
time step and the two solvers agree within 1e-3 of the field peaks, also
while the plasma conducts and charges the dielectric surfaces.

## Limits on the solved variables

Floors and ceilings can be set in `constant/plasmaProperties`. The first
three entries are the defaults for the number densities and the electron
temperature (values shown; the electron temperature has no upper limit unless
`TeMax` is given). A sub-dictionary named after a field sets `min` and/or
`max` for that field alone and takes precedence: `N_<specie>` for the number
densities, `Te`, `Tion` (ion temperature) and `T` (gas temperature).

```
limits
{
    densityFloor    1e4;    // [1/m3]
    TeMin           300;    // [K]
    // TeMax        1e6;    // [K]

    N_electron      { min 1e10; }
    N_Arm           { min 1e12; max 1e20; }
    Te              { min 300; max 1.2e5; }
}
```

The limits are applied after each solve of the field. They are read again
when `plasmaProperties` is edited while the solver runs (with
`runTimeModifiable yes` in `system/controlDict`), so they can be tightened or
relaxed during a run; the log reports the limits in force. Where the electron
density is negligible, for example in a sheath without secondary emission,
the electron energy equation is poorly conditioned; a higher density floor
and a ceiling on `Te` keep the solution bounded there.

## Acceleration of slow neutral species

Metastables reach their periodic steady state over times that are thousands
of periods of the applied voltage. They can be advanced alone, with a large
time step, between blocks of full simulation (`constant/plasmaProperties`):

```
slowSpeciesAcceleration
{
    species         (Arm);
    period          2.5e-8;   // period of the applied voltage [s]
    startTime       1e-5;     // no advances before this time (default 0)
    fullCycles      20;       // periods of full simulation per block
    averageCycles   2;        // periods at the end of a block over which
                              // the chemistry source is averaged
    deltaT          1e-6;     // time step of the advance [s]
    nSteps          20;       // steps per advance
    steadyState     no;       // yes: one steady solve per block instead
    relaxation      1;        // fraction of the change that is applied
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

The advance assumes that the plasma repeats from one period to the next. Use
`startTime` to begin after the discharge has ignited, and keep each advance
(`nSteps*deltaT`) short compared with any slower evolution of the discharge
itself; a long advance or a steady solve otherwise drives the species to a
balance with a plasma state that is still changing.

With `steadyState yes` the averaged equation is solved for its steady state
once per block, and `deltaT` and `nSteps` are optional. Sources that are
nonlinear in the species' own density, such as metastable pooling, are
linearised about the density of the block, so successive blocks act as the
nonlinear iteration; `relaxation` and `maxChangeFactor` keep the changes
moderate. A steady state needs a net loss of the species in every cell; where
there is none, the time-step advance is used if `deltaT` and `nSteps` are
given and the advance is skipped otherwise.

## Spatially varying initial conditions

`setExpressionFields` sets the internal values of scalar fields from formulas
of the cell centre coordinates `x`, `y`, `z` (in metres), given in
`system/setExpressionFieldsDict`:

```
variables               // optional: constants and helper formulas
{
    nBackground 5e14;
    nBlob       2e16;
    x0          1e-3;
    y0          4e-4;
    radius      1.5e-4;
    r2          "sqr(x - x0) + sqr(y - y0)";
}

fields
{
    electron    "nBackground + nBlob*exp(-r2/sqr(radius))";
    Arp1        "nBackground + nBlob*exp(-r2/sqr(radius))";
    Te          "10000 + 5000*x/2e-3";
}
```

Run it in the case directory after the mesh exists and before the solver. Each
listed field must already be in the time directory (default: the start time;
`-time`, `-latestTime` and `-region` are available); its internal values are
replaced and its boundary conditions are kept. The formulas may use
`+ - * / ^`, functions such as `sqr`, `sqrt`, `exp`, `log`, `sin`, `cos`,
`tanh`, `erf`, `min`, `max`, `mag`, `pos`, the constants `pi_` and `e_`, and
the entries of `variables`. `examples/plasma/2DAdaptiveMeshArgon` contains the
dictionary that creates its seeded blob.

## Post-processing

`postprocessing/` holds Python scripts that plot profiles and time histories
of 1D cases (`plot_1d.py`), maps and line cuts of 2D cases including refined
and parallel ones (`plot_2d.py`), and electrode voltage and current
(`plot_electrodes.py`). See `postprocessing/README.md`.

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
