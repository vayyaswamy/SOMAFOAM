# Ambipolar decay of an afterglow (paper case V9)

Electrons at a fixed temperature of 1 eV (`TeMin` and `TeMax` 1 K apart in
`limits`) and argon ions at the gas temperature of 300 K, with constant
mobilities of 30 and 1.44 m2/(V s), Einstein diffusion and no reactions, decay
between two grounded, absorbing walls 2 cm apart. They start in the
fundamental diffusion mode, n0 sin(pi x/L), set by `setExpressionFields`. In
the quasi-neutral limit both densities decay as exp(-t/tau_a), with
tau_a = L^2/(pi^2 D_a) = 28.75 us and the ambipolar diffusion coefficient
D_a = (mu_i D_e + mu_e D_i)/(mu_i + mu_e) = 1.410 m2/s. The walls use the wall
conditions of the plasma examples without secondary emission.

`./Allrun` runs peak densities of 1e17, 1e18 and 1e19 m^-3 on 200 cells and
1e19 m^-3 on 400 and 800 cells, for 55 us each, and `verify.py` compares the
decay rate of the ion inventory over the second half of each run and the final
profile with the ambipolar solution. The five runs took 5.5 min on one core.

Results of 2026-11-17 on Pinnacles:

| n0 [m^-3] | cells | Debye length / L | error of the decay rate | profile error |
|---|---|---|---|---|
| 1e17 | 200 | 2.3e-3 | +8.2 % | 5.0 % |
| 1e18 | 200 | 7.3e-4 | +3.4 % | 2.6 % |
| 1e19 | 200 | 2.3e-4 | +2.1 % | 2.0 % |
| 1e19 | 400 | 2.3e-4 | +1.4 % | 1.3 % |
| 1e19 | 800 | 2.3e-4 | +1.1 % | 0.9 % |

The computed plasma decayed faster than the quasi-neutral solution because the
sheaths at the walls shortened the quasi-neutral region. The error fell with
the Debye length, and at the highest density, where the sheath was thinner
than a cell, with the cell size. Richardson extrapolation of the last three
runs gave about 0.85 % on a fine mesh, which corresponds to sheaths of about
ten Debye lengths at each wall.

The case is a regression test of the bulk transport, the Einstein relation,
the semi-implicit Poisson equation and the wall conditions of the
drift-diffusion model.
