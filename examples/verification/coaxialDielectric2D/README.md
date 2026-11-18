# Coaxial capacitor with a dielectric layer (paper case V14)

The coupled Poisson equation of a plasma region and a dielectric region
across a curved interface, against the exact solution. A quarter of a coaxial
capacitor in the x-y plane has the electrode at r = 1 mm at 10 V, a
dielectric of relative permittivity 5 (region `dielectric_1`) to the
interface at r = 2 mm, and gas to the grounded wall at r = 5 mm, with
symmetry planes on the cuts. The gas holds only the density floor, so there
is no space charge, and the exact potential is logarithmic in each region,
with the charge per length q = 2 pi V0/(ln(b/a)/eps_d + ln(c/b)/eps)
(`verify.py`). `makeCase.py` writes the two block meshes with arc edges,
which match on the interface, and the boundary conditions, starting from
`plasmaDielectric/ArgonDBD`; `dielectricCentres` writes the cell centres of
the dielectric region.

`./Allrun` runs 10, 20 and 40 cells along the arc (with 5, 10 and 20 radial
cells in the dielectric and 15, 30 and 60 in the gas) for five steps of
0.1 ns. The three runs took 16 s on one core.

Results of 2026-11-18 on Pinnacles:

| cells along the arc | max error / V0 |
|---|---|
| 10 | 4.51e-3 |
| 20 | 1.14e-3 |
| 40 | 2.86e-4 |

The potential converged at second order (observed orders 1.99 and 1.99)
across the interface, where the dielectric displacement is continuous and the
permittivity jumps by a factor of five.
