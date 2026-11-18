# One-dimensional adaptive refinement following a moving pulse (paper case V2)

The electron avalanche of `townsendAvalanche1D`, with an initial width of
0.4 mm, on three meshes: uniform with 125 cells (80 um), uniform with 1000
cells (10 um), and 125 base cells with up to three levels of refinement by
`dynamicRefine1DFvMesh` (`constant/dynamicMeshDict.amr`): the relative
gradient of the electron density above 0.02 refines, below 0.005 merges, two
buffer layers, an update every 10 steps. The pulse drifts 3 mm, spreads and
grows over 100 ns with a step of 25 ps, and `verify.py` compares the rates of
growth, drift and spreading of its moments and its final profile with the
exact Gaussian solution. The three runs took 2.5 min on one core.

Results of 2026-11-18 on Pinnacles:

| mesh | cells at the end | growth | drift | diffusion | profile error |
|---|---|---|---|---|---|
| uniform, 125 | 125 | +1.2e-4 | +6.0e-4 | +4.3e-3 | 5.3e-3 |
| uniform, 1000 | 1000 | +1.2e-4 | +1.8e-4 | +3.9e-3 | 1.5e-3 |
| adaptive, 125 base | 630 | +7.4e-5 | +8.4e-5 | +3.3e-3 | 1.9e-3 |

The refinement followed the pulse over 157 updates and kept the finest
cells on its flanks, where the relative gradient was largest. With 63 % of
the cells of the fine uniform mesh, the adaptive run reached nearly its
accuracy, and the excess of the diffusion coefficient in all runs was the
numerical diffusion v^2 dt/2 of the implicit step (3.75e-3 of D).
