# Argon global model at set power (paper case B1)

The zero-dimensional model (`plasmaModel zeroD`) of an argon discharge in a
cylinder of radius 0.15 m and length 0.3 m, driven at a set absorbed power
(`zeroDDrive`, mode `power`, with the coupled implicit step) and with the wall
losses of `zeroDWall`, against the global model of Lieberman and Lichtenberg
(Principles of Plasma Discharges and Materials Processing, 2nd ed., 2005,
chapter 10). The reactions are elastic collisions, ionization and excitation
of the ground state with the Maxwellian rate coefficients of that model,
which `makeRates.py` writes as tables; the excitation is an energy loss only,
so that there are no metastables. In the global model the particle balance
N K_iz(Te) = u_B A_eff/V gives Te, and the power balance of the electrons
with the collisional losses and 2 kB Te per electron lost gives n_e
(`compare.py`, which evaluates the fits directly).

`./Allrun` runs 100, 300 and 1000 W at 10 mTorr and 3, 10 and 30 mTorr at
300 W, each for 10 ms with a step of 1 us, which reaches the steady state.

Results of 2026-11-18 on Pinnacles:

| P (W) | p (mTorr) | Te (eV) | rel. error | n_e (m^-3) | rel. error |
|---|---|---|---|---|---|
| 100 | 10 | 2.2541 | -4.2e-5 | 5.596e16 | +5.7e-4 |
| 300 | 10 | 2.2542 | +2.9e-6 | 1.679e17 | +6.3e-4 |
| 1000 | 10 | 2.2545 | +1.6e-4 | 5.597e17 | +8.4e-4 |
| 300 | 3 | 2.8930 | +1.3e-4 | 1.272e17 | +7.2e-4 |
| 300 | 30 | 1.8664 | -4.5e-5 | 2.095e17 | +5.8e-4 |

The electron temperature agreed to 2e-4 and the electron density to 1e-3,
with n_e higher than the global model by 6e-4 to 8e-4 in all runs; the
origin of this small common offset has not been traced. The electron
density was proportional to the power and the temperature independent of
it, and the temperature fell with the pressure, as in the global model.
