# SOMAFOAM user guide

This guide explains how to build SOMAFOAM, run its examples, run the test suite and
submit any of these as a batch job on a SLURM cluster. Every task goes through one
command, `scripts/somafoam`, which loads the compiler modules and the SOMAFOAM
environment by itself, so nothing has to be sourced beforehand. The equations, schemes,
boundary conditions and dictionary entries are described in the reference manual
[`doc/manual/somafoam_manual.pdf`](manual/somafoam_manual.pdf), and the ion distribution
workflow in [`doc/manual/ion_distribution_workflow.pdf`](manual/ion_distribution_workflow.pdf).

## Contents

1. [Quick start](#1-quick-start)
2. [Building](#2-building)
3. [Running an example](#3-running-an-example)
4. [Running the tests](#4-running-the-tests)
5. [Running on a SLURM cluster](#5-running-on-a-slurm-cluster)
6. [Post-processing](#6-post-processing)
7. [Setting up a new case](#7-setting-up-a-new-case)
8. [Adding a test case](#8-adding-a-test-case)
9. [Where the scripts are](#9-where-the-scripts-are)
10. [Troubleshooting](#10-troubleshooting)

## 1. Quick start

```
git clone https://github.com/vayyaswamy/SOMAFOAM
cd SOMAFOAM
scripts/somafoam build                       # about 30 min on 28 cores
scripts/somafoam test flow                   # three flow tests, about 1.5 min
scripts/somafoam list                        # all examples
scripts/somafoam run zeroD/ArgonGlobal       # copy an example to ./runs and run it
```

To type `somafoam` instead of `scripts/somafoam`, add the scripts directory to the
search path in `~/.bashrc`.

```
export PATH=$HOME/SOMAFOAM/scripts:$PATH     # adjust to where the clone is
```

On a cluster, put `submit` in front of `build`, `run` or `test` to run the same command
as a batch job, as described in Section 5.

```
somafoam submit build
somafoam submit test
somafoam submit run plasma/2DAdaptiveMeshArgon --np 8
```

## 2. Building

SOMAFOAM needs gcc and g++, make, flex, Python 3 and OpenMPI with its development headers.
It has been built with gcc 10.2.1 and OpenMPI 4.1.0 on Debian 11 and with gcc 11.5.0 and
OpenMPI 4.1.8 on the Pinnacles cluster at UC Merced. The repository holds sources only,
and the libraries and applications have to be built on each machine because they depend
on its compiler and C library.

```
somafoam build            # all cores
somafoam build -j 8       # 8 compile processes
```

The command builds wmake, the third-party libraries, all libraries in `src/` and all
applications in `applications/`, and writes the output to `logs/build-<date>.log`. It
then checks that `somaFoam`, `blockMesh` and `libplasmaCookBook.so` were installed and
reports success or the first error lines of the log. The libraries go to
`lib/linux64GccDPOpt` and the applications to `bin/linux64GccDPOpt`. A full build takes
about 30 minutes on 28 cores, and a rebuild after a change compiles only what changed.

To use the applications directly in a shell, without the `somafoam` command, load the
environment first.

```
source scripts/sites/pinnacles.sh     # only on a cluster that needs modules
source etc/bashrc
somaFoam -help
```

A build made on one machine does not run on a machine with an older C library, so a
build from a cluster usually fails on a workstation with a message about `GLIBC_2.34`.
Build on each machine where you run.

## 3. Running an example

The examples are in `examples/`, grouped by the kind of problem. `somafoam list` prints
every example with its solver and whether it is a test.

| Directory | Content |
|---|---|
| `flow/` | gas flow with analytic solutions, which are tests |
| `plasma/` | one- and two-dimensional discharges, adaptive refinement, external circuits, ion distributions |
| `plasmaDielectric/`, `dielectricCases/` | discharges with dielectric regions |
| `zeroD/` | global models and chemistry on one cell |
| `electronBoltzmann/` | electron transport and rate tables from cross sections |

`somafoam run` copies an example into a run directory and runs it there, so the example
in the repository stays clean.

```
somafoam run plasma/2DAdaptiveMeshArgon                   # serial
somafoam run plasma/2DAdaptiveMeshArgon --np 8            # parallel on 8 processes
somafoam run plasma/2DAdaptiveMeshArgon --end-time 1e-9   # a few steps to check the setup
somafoam run 2DAdaptive --dry-run                         # show the steps only
somafoam run zeroD/ArgonGlobal --tag pulse2               # second copy, runs/zeroD_ArgonGlobal_pulse2
```

The case can be named by its path under `examples/`, by any part of the path that is
unique, or by the path of any case directory on disk. The run directory is
`./runs/<group>_<case>` by default. `--dir` or the variable `SOMAFOAM_RUN_DIR` sets another
parent directory, which on a cluster should be on the scratch file system.

A case with an `Allrun` script is run by that script. Any other case is run with the
standard steps, which are `Allmesh` or `prepareCase` if the case has one, otherwise
`blockMesh` when the mesh is not included, then `setExpressionFields` when the case has a
`setExpressionFieldsDict`, and finally the solver named in `system/controlDict`. With
`--np N` the copy is decomposed into N parts with `decomposePar`, including each dielectric
region, and the solver runs with `mpirun -np N ... -parallel`. Every step writes its output
to `log.<step>` in the run directory, and when a step fails the end of its log is printed.
`--dry-run` lists the steps without running them, and `--in-place` runs in the given
directory instead of a copy.

`somafoam clean DIR` removes the time directories, processor directories and logs of a
run and keeps `0/`, `constant/` and `system/`.

## 4. Running the tests

An example with an `Allrun` script is a test. Its `Allrun` builds the mesh, runs the
solver and compares the result with a reference, which is an analytic solution, results of
another code or published data. The comparison prints the error and ends with PASS or
FAIL.

```
somafoam test                       # all tests
somafoam test flow icp              # tests whose path contains flow or icp
somafoam test --list                # names only
somafoam test --short               # shortened runs, if the Allrun supports them
somafoam test --jobs 4 --keep       # four at a time, keep the run directories
somafoam test --scratch /scratch/$USER/tests
```

Each test runs in a copy under the scratch directory, which is `--scratch`,
`SOMAFOAM_SCRATCH` or the system temporary directory, in that order. On a cluster the
temporary directory is local to the compute node and disappears with the job, so use
`--scratch` together with `--keep` when the run directories are needed afterwards. The
summary lists each test as PASS, FAIL, RAN or ERROR. RAN means that the case finished but
has no reference yet, and ERROR means that a step failed or the test timed out after
`--timeout` seconds, 7200 by default. The full output of each test is in `log.Allrun` in
its run directory, and the summary is also written to `test/results.json`. The command
returns a non-zero exit code when any test failed, so it can be used in scripts.

The cases used in the SOMAFOAM 2.0 paper, and the example directory of each, are listed
in [`examples/README.md`](../examples/README.md).

## 5. Running on a SLURM cluster

`somafoam submit` writes a batch script for `build`, `run` or `test` and submits it with
`sbatch`. Everything after the subcommand is passed on unchanged.

```
somafoam submit build
somafoam submit run plasma/2DCurvedElectrodeReactor --np 28 --dir /scratch/$USER/runs
somafoam submit test --jobs 4 --scratch /scratch/$USER/tests --keep
somafoam submit run zeroD/CF4O2 --partition test --time 00:30:00
somafoam submit run plasma/100TorrArgonPlasma --sbatch=--mem=64G
somafoam submit build --print       # show the job script without submitting
```

| Option | Meaning |
|---|---|
| `--partition P` | SLURM partition, by default the one in the site file |
| `--time HH:MM:SS` | wall-time limit, by default the one in the site file |
| `--cpus N` | cores for `build` and `test`; `run --np N` asks for N tasks instead |
| `--sbatch=...` | any further `#SBATCH` line, repeatable |
| `--print` | print the job script and do not submit it |

The job script and the job output go to `logs/` in the repository, and the output file is
named `sf-<command>-<case>-<jobid>.out`. The job starts in the directory from which it was
submitted, so a relative `--dir` or the default `./runs` refers to that directory. Follow
the jobs with `squeue -u $USER` and cancel one with `scancel <jobid>`.

### Site files

A site file in `scripts/sites/` loads the modules a machine needs and sets the defaults of
`submit`. The site is taken from the variable `SOMAFOAM_SITE`, or else recognised from the
SLURM cluster name or the host name, and `scripts/sites/default.sh` is used when nothing
matches. `somafoam env` shows the site in use, the solver and MPI found, and the default
directories.

```
# scripts/sites/pinnacles.sh
module purge
module load openmpi/4.1.8-gcc-11.5.0
SLURM_PARTITION=short          # test 1 h, short 6 h, medium 1 day, long 3 days
SLURM_TIME=06:00:00
SLURM_CPUS=28
```

A new cluster needs a copy of this file named after the cluster, with its own modules and
partitions. The same compiler and MPI have to be loaded for the build and for the runs,
which the site file guarantees.

### Notes for Pinnacles

The `/data` file system has a quota of 500 GB with a hard limit of 512 GB per user, and
`/scratch` is large but purged. Builds and source trees fit on `/data`, whereas runs and
test directories belong on `/scratch/$USER`. Each login node has its own `/tmp`, so a file
written to `/tmp` in one session may not be visible in the next. In a large source tree
`git status -uno` is much faster than `git status`, which scans every untracked file.

## 6. Post-processing

The scripts in `postprocessing/` read SOMAFOAM cases directly and plot profiles, time
histories, two-dimensional maps and electrode currents without ParaView. They need
`numpy` and `matplotlib`, and the case must be written in ASCII. Their use is described in
[`postprocessing/README.md`](../postprocessing/README.md).

```
postprocessing/plot_1d.py runs/plasma_100TorrArgonPlasma --fields electron Arp1 Te --times latest
postprocessing/plot_2d.py runs/plasma_2DAdaptiveMeshArgon --fields electron Te --time latest -o ne.png
postprocessing/plot_electrodes.py runs/plasma_2DAdaptiveMeshArgon
```

ParaView reads a case through an empty file named `<case>.foam` in the case directory, and
`foamToVTK` converts a case for other viewers. Parallel runs are combined with
`reconstructPar`, or with the `reconstructParScript` of cases with dielectric regions.

## 7. Setting up a new case

The simplest start is a copy of the example closest to the new problem, made with
`somafoam run EXAMPLE --dry-run` to see its steps and then with
`cp -r examples/<group>/<case> mycase`. A case has three directories. The `0/` directory
holds the initial and boundary fields, with one file per species, `Te`, `Phi` and, for
flow, `U`, `p` and `T`. The `constant/` directory holds the mesh, `plasmaProperties` with
the plasma model and species, `chemistryProperties` and `reactions` with the chemistry,
the transport tables, `electroMagnetics` with the fields and electrodes and
`regionProperties` with the plasma, dielectric and electromagnetic regions. The `system/`
directory holds `controlDict` with the time step and output, `fvSchemes`, `fvSolution`,
`blockMeshDict` and `decomposeParDict`. Every entry is explained in the reference manual.
After editing, run the case in place.

```
somafoam run mycase --in-place
somafoam submit run mycase --in-place --np 16
```

## 8. Adding a test case

A test is an example directory with an `Allrun` script, a script that compares the result
with the reference and prints PASS or FAIL as its last line, a `reference/` directory with
the reference data and its source, and a `README.md` that states what the case tests, where
the reference comes from and the paper section and figure it belongs to. The case
`examples/flow/channelStartupProfile` is a short complete example.

```
#!/bin/sh
cd ${0%/*} || exit 1
blockMesh > log.blockMesh 2>&1 || exit 1
somaFoam  > log.somaFoam  2>&1 || exit 1
python3 ../verify.py channelStartupProfile .
```

An `Allrun` may shorten the run when `SOMAFOAM_TEST_SHORT=1`, which `--short` sets, and
must exit with a non-zero code when a step fails so that the test is reported as an ERROR
rather than a FAIL. A new test is found by `somafoam test` without any registration.

## 9. Where the scripts are

| Path | Purpose |
|---|---|
| `scripts/somafoam` | the command described in this guide |
| `scripts/sites/*.sh` | modules and SLURM defaults of each machine |
| `install.sh` | the build itself, called by `somafoam build` |
| `test/runExamples.py` | the test driver, called by `somafoam test` |
| `examples/*/*/Allrun`, `verify.py`, `compare.py` | the run and comparison of each test |
| `examples/*/*/prepareCase`, `caseClean` | mesh and region set-up, and clean-up, of examples with dielectrics |
| `postprocessing/*.py` | plots read directly from the case files |
| `etc/bashrc` | the SOMAFOAM environment |
| `logs/` | build logs, job scripts and job output, not tracked by git |
| `runs/` | default run directories, not tracked by git |

## 10. Troubleshooting

| Message or symptom | Cause and remedy |
|---|---|
| `somaFoam: command not found` in a shell | the environment is not loaded; use `somafoam run`, or source the site file and `etc/bashrc` |
| `GLIBC_2.34 not found` | the build comes from a newer system; build on this machine |
| `Disk quota exceeded` | the file system is full; move runs to scratch and remove old time directories with `somafoam clean` |
| a test reports ERROR | read `log.Allrun` and the step logs in the kept run directory, `somafoam test NAME --keep --scratch DIR` |
| a parallel run stops at start-up | `numberOfSubdomains` must match `--np`; `somafoam run --np` sets it in the copy |
| `run` refuses because the run directory exists | add `--force` to replace it or `--tag NAME` for another copy |
