# Electron energy relaxation in a global model (paper cases V13 and V7)

A zero-dimensional argon gas at 1 Torr and 300 K (`plasmaModel zeroD`) with a
constant electron mobility of 30 m2/(V s), no reactions and a reduced field of
2.75 Td switched on for 20 us and off for 30 us (`zeroDDrive`). The electron
energy equation is then

    d/dt (3/2 kB Te) = e mu E^2 - 3 kB (m/M) nu (Te - T),   nu = e/(m mu),

so Te rises towards T + mu^2 E^2 M/(3 kB) = T + 11294 K at the rate
2 e/(M mu) = 1/(6.21 us) and decays back to T after the field is switched
off. With `writeEnergyBudget yes` the terms of the energy equation are written
at the write times, and their residual checks the closure of the budget.

`./Allrun` runs time steps of 10 and 40 ns, and `verify.py` compares Te with
the exact solution and checks the budget. The two runs took 80 s on one core.

Results of 2026-11-18 on Pinnacles:

| time step | max error of Te / rise | budget residual / heating |
|---|---|---|
| 10 ns | 1.60e-3 | 1.2e-13 |
| 40 ns | 6.41e-3 | 1.4e-14 |

The error fell at first order with the time step (observed order 1.00), as
expected of the implicit Euler step, and the budget closed to round-off.

Note: the reaction list must not contain electron-impact elastic reactions
when their rates are not available. With such a reaction in the list and
`chemistry off`, the elastic loss takes its collision frequency from that
reaction, whose rate is then zero, and the electrons lose no energy to the
gas. The first set-up of this case did so and heated without bound.
