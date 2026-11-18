# Electron avalanche in a uniform field (paper case V12)

A Gaussian pulse of electrons and ions (peak 1e10 m^-3, width 0.2 mm, at
x = 7 mm) in argon at 1 Torr between electrodes 1 cm apart at 10 V, with a
constant electron mobility of 30 m2/(V s), the electron temperature held at
0.1 eV and a single ionization reaction of constant rate (Arrhenius A with no
temperature dependence) that gives nu = 1.000e7 1/s. The space charge of the
pulse is negligible, so the field is uniform, and the electron density obeys
dn/dt + d(v n)/dx - D d2n/dx2 = nu n with v = -mu E = -3e4 m/s and
D = mu kB Te/e = 3.001 m2/s. The exact solution is a Gaussian whose content
grows as exp(nu t), whose centre moves at v and whose variance grows as 2 D t.

`./Allrun` runs time steps of 100 and 25 ps over 100 ns on 1000 cells, and
`verify.py` compares the rates of growth, drift and spreading of the moments
of the computed density with the exact ones. The two runs took 1 min 46 s on
one core.

Results of 2026-11-18 on Pinnacles:

| time step | ionization frequency | drift velocity | diffusion coefficient | v^2 dt/(2D) |
|---|---|---|---|---|
| 100 ps | +5.0e-4 | +9.6e-4 | +1.58e-2 | 1.50e-2 |
| 25 ps | +1.3e-4 | +2.1e-4 | +3.81e-3 | 3.75e-3 |

The excess of the diffusion coefficient was the numerical diffusion v^2 dt/2
of the implicit Euler step and fell at first order with the step (observed
order 1.03).

The pressure in `0/p` must be consistent with the gas density in `0/Ar`: the
chemistry forms the concentrations from the gas density of the
thermophysical model, which follows from the pressure. The first set-up of
this case took the pressure of 100 Torr from `plasma/ExternalCircuitDC`, whose
`0/Ar` holds the density of 1 Torr, and ionized 100 times too fast.
