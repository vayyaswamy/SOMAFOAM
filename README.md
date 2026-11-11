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

## Features

Options are set in the case dictionaries; the examples show them in use.

| Feature | Where |
|---|---|
| Plasma with dielectric regions: `somaFoam` (coupled) or `plasmaMultiRegionFoam` (iterated) | `examples/plasmaDielectric/ArgonDBD` |
| Scharfetter-Gummel fluxes: `fluxScheme scharfetterGummel;` in `constant/plasmaProperties` | `examples/plasma/1TorrArgonPlasma` |
| Inner iterations, limits on the variables, acceleration of slow species: `innerIterations`, `limits`, `slowSpeciesAcceleration` in `constant/plasmaProperties` | `examples/plasma/100TorrArgonPlasma` (acceleration only) |
| Adaptive mesh refinement (1D, 2D, 3D): `constant/dynamicMeshDict`; `rebalancePar` for parallel runs | `examples/plasma/2DAdaptiveMeshArgon` |
| Initial fields from formulas of x, y, z: `setExpressionFields` | `examples/plasma/2DAdaptiveMeshArgon/system/setExpressionFieldsDict` |
| Stop at a steady or periodic state when window-averaged fields stop changing: `convergenceStop` in `system/controlDict` (fields, window, tolerance, nWindows, startTime) | `applications/somaFoam/convergenceStop.H` |
| Electrode voltage and current: `electrodeVoltageCurrent` function object | `system/controlDict` of the examples |
| RF electrode at a given mean power (constant or a profile in time) with a blocking capacitor and sine, multi-harmonic or tabulated voltage waveforms: `rfPowerElectrode` condition for `Phi` | `doc/manual` |
| Electron transport and rate tables from cross sections (LXCat format), two-term or Monte Carlo: `electronBoltzmann`; run by `somaFoam` at start-up with `electronBoltzmann yes;` in `constant/plasmaProperties` | step-by-step guide: [`applications/preProcessing/electronBoltzmann/README.md`](applications/preProcessing/electronBoltzmann/README.md); case `examples/electronBoltzmann/Argon` |
| Ion energy and angular distributions at the walls from the fields of a run (Monte Carlo ion tracking): `ionTracker` | `examples/plasma/HeliumIonDistribution`, `doc/manual/ion_distribution_workflow.pdf` |
| Gas flow driven by the plasma, with the equations of `rhoSimpleFoam` (steady, at the start) or `rhoPimpleFoam` (transient) and the compressible turbulence models; species carried by the gas: `constant/flowProperties` | `examples/flow` (verification cases with analytic solutions and the known limitations), `applications/somaFoam/flow/createFlowFields.H` |
| Zero-dimensional (global) model on one cell driven by a prescribed E/N(t) or Te(t): `plasmaModel zeroD<constGasThermoPhysics>;` with `zeroDDrive`, `plasmaChemistryModel temporal;` and `poissonEquationSolver none;` | `examples/zeroD/ArgonGlobal` |
| Plasma chemistry evaluated field by field (faster; same reactions file): `psiChemistryModel plasmaFieldChemistryModel<constGasThermoPhysics>;` in `constant/chemistryProperties` | `src/thermophysicalModels/chemistryModel/chemistryModel/plasmaFieldChemistryModel` |
| Electron energy as the energy density 3/2 n_e kB Te with the energy flux in drift-diffusion form (Hagelaar), conservative by construction: `electron eEnergyDensity;` in `energyModel` of `constant/plasmaProperties` (default `efullImplicit`, Te) | `src/plasmaCookBook/plasmaEnergyModels/eTemperature/eEnergyDensity` |
| Drift fluxes through the faces from the potential difference across each face, `driftFlux face;` in `constant/plasmaProperties` (default `cell`, interpolated drift velocity) | `src/plasmaCookBook/multiSpeciesPlasmaModels/multiSpeciesPlasmaModel/multiSpeciesPlasmaModelNumerics.C` |
| Electron wall conditions (`driftDiffusionElectronDensity`, `electronTemperature`): thermal speed of the wall flux from the wall-cell temperature, `thermalSpeed cell;` (default `wall`), and ion flux of secondary emission computed at the wall from the ion densities, `emissionIonFlux wall; emissionIons (Arp1);` (default `field`) | `src/plasmaCookBook/plasmaBoundaryConditions/plasmaWallFluxes.H` |
| Interpolation of tabulated rates: `rateInterpolation logarithmic;` (default, ln k linear in Te) or `linear` (k linear in Te) in `constant/chemistryProperties` | `src/thermophysicalModels/specie/reaction/reactionRate/TabularReactionRate` |
| Faster time steps: reaction rates of the species equations `speciesChemistry perSpecies;` (default), `perPass` or `perStep` (rates of the start of the step), and `electronTemperature 1;` (one pass), in `innerIterations` of `constant/plasmaProperties`. With `perStep`, one pass, `Gauss linear` convection and `runTimeModifiable no;`, a 1D argon CCP takes as long per step as a single-pass time-marching solver | `src/plasmaCookBook/multiSpeciesPlasmaModels/multiSpeciesPlasmaModel/multiSpeciesPlasmaModelNumerics.C` |
| Python plotting scripts | `postprocessing/README.md` |

Contributors:
1) Venkattraman Ayyaswamy (https://me.ucmerced.edu/content/venkattraman-venkatt-ayyaswamy)
2) Abhishek Kumar Verma 
3) Saurav Gautam (https://me.ucmerced.edu/content/saurav-gautam)
4) Jose Alfredo Millan Higuera (https://www.linkedin.com/in/jmillanhiguera/)

For research conducted using SOMAFOAM please cite this paper: https://doi.org/10.1016/j.cpc.2021.107855

Paper that compares SOMAFOAM results with experiment: 
https://doi.org/10.1063/5.0041386
