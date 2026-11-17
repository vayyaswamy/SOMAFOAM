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
| Manufactured solution, order of accuracy | `verification/manufacturedSolution2D` (planned) | exact |
| Moving ionization front, adaptive refinement | `verification/ionizationFront1D` (planned) | traveling wave |
| Gas flow: advection-diffusion, channel start-up, transient Poiseuille profile | `flow/advectionDiffusion`, `flow/channelStartup`, `flow/channelStartupProfile` | analytic |
| Flowing afterglow | `flow/flowingAfterglow` (planned) | eigenfunction series |
| Inductive field in a uniform plasma | `icp/skinEffectCylinder`, `icp/skinEffectSlab` (planned) | Bessel and skin-depth solutions |
| Serial and parallel adaptive runs | `plasma/2DAdaptiveMeshArgon` | identical results |
| Argon global model, steady and pulsed | `zeroD/ArgonGlobal`, `zeroD/ArgonPulsed` (planned) | published global models |
| CF4/O2 chemistry | `zeroD/CF4O2` | Cantera |
| Electron transport tables | `electronBoltzmann/Argon` | BOLSIG+ |
| Low-pressure CCP benchmark | `plasma/TurnerBenchmark` (planned) | Turner et al. 2013 |
| Positive streamer in air | `streamer/BagheriCase1` to `Case3` (planned) | Bagheri et al. 2018 |
| Helium DBD and period multiplication | `dielectricCases/sinusoidal-helium`, `plasmaDielectric/heliumPeriodMultiplication` (planned) | published fluid results |
| Ion distributions | `plasma/HeliumIonDistribution` | WarpX PIC |
| Planar-coil argon ICP | `icp/GodyakPlanarCoil` (planned) | Godyak et al. 2002 |
| GEC ICP reference cell | `icp/GECICPCell` (planned) | Miller et al. 1995 |
| Surface DBD actuator | `flow/surfaceDBDActuator` (planned) | published power, force, and velocity |
| Positive streamer against experiment | `streamer/airExperiment` (planned) | published measurements |

The other examples (`plasma/100TorrArgonPlasma`, `plasma/1TorrArgonPlasma`,
`plasma/2DCurvedElectrodeReactor`, `plasma/ExternalCircuitDC`, `plasmaDielectric/ArgonDBD`)
show features and are not tests.
