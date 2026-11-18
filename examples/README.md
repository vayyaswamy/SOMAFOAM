# Examples

Every example is a complete case. Examples with an `Allrun` script are also tests: `Allrun`
makes the mesh, runs the solver, and compares the result with its reference, printing the
error and PASS or FAIL. Run all of them, or a subset, with

    source etc/bashrc
    python3 test/runExamples.py            # all examples with an Allrun
    python3 test/runExamples.py flow icp   # paths containing flow or icp
    python3 test/runExamples.py --short    # short versions for quick checks

`--short` sets `SOMAFOAM_TEST_SHORT=1`, which an `Allrun` may use to run fewer steps.

## Layout of a test case

| File | Content |
|---|---|
| `Allrun` | mesh, run, post-processing, comparison |
| `verify.py` or `compare.py` | error against the reference, last line PASS or FAIL |
| `reference/` | analytic solution, benchmark data, or digitized published data with its source |
| `README.md` | what the case tests, the reference, and the paper section and figure |

## Cases of the SOMAFOAM 2.0 paper

Release tag `v2.0.0`. Entries marked *planned* are being added for the release.

| Paper case | Directory | Reference |
|---|---|---|
| Manufactured solution, order of accuracy | `verification/manufacturedSolution2D` | manufactured solution; second order on a uniform mesh, first order on a distorted one |
| Adaptive refinement following a moving pulse | `verification/adaptiveFront1D` | Gaussian solution; 630 adaptive cells nearly as accurate as 1000 uniform |
| Electron energy relaxation and energy budget | `verification/electronEnergyRelaxation0D` | exact relaxation; 1.6e-3 of the rise, first order in the step, budget to 1e-13 |
| Electron avalanche in a uniform field | `verification/townsendAvalanche1D` | Gaussian solution; growth and drift to 2e-4, numerical diffusion v^2 dt/2 |
| Coaxial capacitor with a dielectric layer, curved interface | `verification/coaxialDielectric2D` | logarithmic solution; second order, 2.9e-4 of the voltage |
| Ambipolar decay of an afterglow | `verification/ambipolarDecay1D` | ambipolar diffusion; error falls with the Debye length and the cell size, 1.1 % on the finest run |
| External circuit, RC charging of an empty gap | `verification/externalCircuitRC` | backward Euler and exact RC solutions; 2e-5 of the source |
| Gas flow: advection-diffusion, channel start-up, transient Poiseuille profile | `flow/advectionDiffusion`, `flow/channelStartup`, `flow/channelStartupProfile` | analytic |
| Flowing afterglow | `flow/flowingAfterglow` (planned) | eigenfunction series |
| Inductive field in a uniform plasma | `verification/icpSkinEffect` | Bessel solution; second order, power to 3e-5 on a 0.25 degree wedge |
| Serial and parallel adaptive runs | `plasma/2DAdaptiveMeshArgon` | serial run; round-off on a fixed mesh, 1e-4 in current with refinement |
| Argon global model, steady and pulsed | `zeroD/ArgonGlobal`, `zeroD/ArgonPulsed` (planned) | published global models |
| CF4/O2 chemistry | `zeroD/CF4O2` | Cantera |
| Electron transport tables | `verification/boltzmannModelGas` | exact Davydov and Maxwellian distributions of model gases; order 1.5, 5e-4 on 800 cells |
| Low-pressure CCP benchmark | `plasma/TurnerBenchmark` (planned) | Turner et al. 2013 |
| Positive streamer in air | `streamer/BagheriCase1` to `Case3` (planned) | Bagheri et al. 2018 |
| Helium DBD and period multiplication | `dielectricCases/sinusoidal-helium`, `plasmaDielectric/heliumPeriodMultiplication` (planned) | published fluid results |
| Ion distributions | `plasma/HeliumIonDistribution` | WarpX PIC |
| Planar-coil argon ICP | `icp/GodyakPlanarCoil` (planned) | Godyak et al. 2002 |
| GEC ICP reference cell | `icp/GECICPCell` (planned) | Miller et al. 1995 |
| Surface DBD actuator | `flow/surfaceDBDActuator` (planned) | published power, force, and velocity |
| Positive streamer against experiment | `streamer/airExperiment` (planned) | published measurements |

The feature examples (`plasma/100TorrArgonPlasma`, `plasma/1TorrArgonPlasma/defaultSchemes`
and `scharfetterGummel`, `plasma/2DCurvedElectrodeReactor`, `plasma/ExternalCircuitDC`,
`plasmaDielectric/ArgonDBD`) are regression tests: their `Allrun` runs the first 50 to 200
steps in `runs/r` and compares the potential, the electron density and Te with the values
stored in `reference/fields.json` to a relative 1e-6 (`test/regression.py`); `./Allrun store`
writes new references after an intended change of the results. The case directories
themselves are not changed, so the examples run as before with `somafoam run`.
