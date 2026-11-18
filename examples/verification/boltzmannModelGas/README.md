# Two-term Boltzmann solver against exact distributions (paper case V15)

`electronBoltzmann` against the exact solutions of the two-term equation for
model gases with elastic collisions only (mass ratio of helium, 300 K). The
zero-flux solution is the Davydov distribution,
F0 = C exp(-int de/(kB T/e + (M/6m) (E/N)^2/(e sigma^2))), which for a
constant cross section of 1e-19 m2 is the Druyvesteyn distribution at T = 0,
and for sigma = 1e-19 (1 eV/eps)^1/2, a constant collision frequency, a
Maxwellian. `exact.py` writes the cross sections (LXCat format, 1e-4 to
1000 eV) and computes the exact mean energy and mu N by quadrature.

`./Allrun` runs 0.5, 1, 2 and 4 Td on 200, 400 and 800 energy cells for
both gases, and `verify.py` prints the relative errors.

Results of 2026-11-18 on Pinnacles, relative errors of the mean energy and of
mu N on 800 cells and the observed order:

| gas | E/N (Td) | mean energy (eV) | error | mu N error | order |
|---|---|---|---|---|---|
| constant sigma | 0.5 | 0.2051 | -1.1e-4 | -1.1e-4 | 1.50 |
| constant sigma | 4 | 1.4809 | -6.0e-5 | -6.0e-5 | 1.50 |
| constant nu | 0.5 | 0.0844 | -4.2e-4 | -3.8e-4 | 1.46 |
| constant nu | 4 | 2.9573 | -4.8e-4 | -4.3e-4 | 1.46 |

The mean energy and the mobility converged to the exact values at an
observed order of 1.5 in the number of energy cells, with errors below 5e-4
on 800 cells.
