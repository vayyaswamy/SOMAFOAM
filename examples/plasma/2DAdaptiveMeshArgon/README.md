# 2D adaptive mesh refinement example

Argon at 100 Torr between a driven electrode (x = 0, 50 V + 100 V sin(2 pi 40 MHz t))
and a grounded electrode (x = 2 mm), on a 2 mm x 0.8 mm two-dimensional mesh
with symmetry planes on the sides. A round blob of plasma is seeded in the
middle of the gap, so the refined region is not aligned with the grid.

- Base mesh 75 x 30 cells, one level of adaptive refinement
  (`constant/dynamicMeshDict`: `dynamicPolyRefinementFvMesh` with the
  `plasmaIndicatorRefinement` selection on the relative gradient of the
  electron and ion densities).
- `driftDiffusion` model, fixed time step 1e-10 s, 0.1 microseconds (1000 steps).
- Non-orthogonal correction is switched on in `system/fvSchemes`
  (`corrected` Laplacian and surface-normal gradient schemes); it roughly
  halves the error at the fine/coarse interfaces.

Run (serial; about 25 minutes):

```
blockMesh      # only needed after caseClean; the mesh is included
somaFoam
```

The mesh in each time directory (`<time>/polyMesh`) is the refined one, and
`refinementIndicator` shows where refinement was requested. Electrode voltage
and current are written to `electrodes/0/<patch>.dat`.

Compared with a uniform 150 x 60 mesh at 0.1 microseconds this case agrees
within 2 % of the peak in electron density (0.4 % RMS).

To try two levels, set `maxRefinementLevel 2;` (the cell count can then reach
36000) and consider `lowerRefineLevel 0.2; unrefineLevel 0.07;`.
