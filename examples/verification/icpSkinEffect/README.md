# Inductive field of a uniform plasma cylinder (paper case V5)

The field solve of the inductive (ICP) option of `somaFoam` against the exact
solution of a uniform plasma cylinder inside a long solenoid. The test
application `icpFieldTest` (`applications/verification/icpFieldTest`) builds
the coil current and solves for the azimuthal vector potential with the
functions of `applications/somaFoam/inductiveField.H`, which `somaFoam` uses,
with the cold-plasma conductivity of a uniform electron density instead of a
plasma solution.

The plasma (n0 = 1e17 m^-3, nu_m = 1e7 1/s, 13.56 MHz, skin depth 1.69 cm)
fills r < 5 cm. The solenoid is 20 contiguous square turns of 2.5 mm at
r = 7 cm, 1 A per turn, and a conducting wall at r = 10 cm sets A = 0. The
mesh is an axisymmetric wedge, 50 mm long with zero-gradient ends, so the
solenoid is infinitely long and the field is uniform along z. `makeCase.py`
writes the mesh and the turns (`coilTurns`). With kappa^2 = i w mu0 sigma the
exact axial field is B1 I0(kappa r)/I0(kappa Rp) in the plasma, uniform
between the plasma and the winding, falls linearly across the winding and is
uniform outside it, and B1 follows from the zero flux inside the wall
(`verify.py`).

`./Allrun` runs 160, 320 and 640 radial cells on a wedge of 0.25 degrees, on
which the edges of the plasma and of the winding are faces, and 640 cells on
wedges of 1 and 5 degrees. The whole run took 8 s on one core.

Results of 2026-11-17 on Pinnacles:

| radial cells | wedge | max error of A | error of the power |
|---|---|---|---|
| 160 | 0.25 deg | 1.29e-3 | 3.3e-4 |
| 320 | 0.25 deg | 3.19e-4 | 8.8e-5 |
| 640 | 0.25 deg | 7.66e-5 | 2.8e-5 |
| 640 | 1 deg | 7.60e-5 | 1.5e-4 |
| 640 | 5 deg | 1.73e-3 | 3.3e-3 |

The vector potential converged at second order (observed orders 2.01 and
2.06), and the field was uniform along z to 1e-12. The flat faces of a wedge
misrepresent the cell volumes and the radius of the cell centres by amounts
of the order of the square of the wedge angle, which gave errors of 3.3e-3 in
the power on the 5 degree wedge and 1.5e-4 on the 1 degree wedge. ICP cases
should therefore use wedges of 1 degree or less when the absorbed power is
needed to better than about 0.1 %.
