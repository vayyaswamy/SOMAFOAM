# Gas flow: verification cases

Cases with analytic solutions that check the gas flow of `somaFoam`
(`constant/flowProperties`) and its coupling to the species. The flow is solved
with the equations of the foam-extend solvers `rhoSimpleFoam` (steady) and
`rhoPimpleFoam` (transient) and the compressible turbulence models
(`constant/turbulenceProperties`, laminar here). Their dictionaries are in
`system/fvSolution`: `SIMPLE` (with `rhoMin`, `rhoMax`, `pRefCell`,
`pRefValue`) or `flowPIMPLE` (with `pMin`), the relaxation factors, and the
solvers of `U`, `p`, `hs` and `rho`.

There is no discharge in the cases (no voltage, `chemistry off`, a trace of
charged species), so only the flow and the transport of the species by the gas
are tested. The gas is argon at 1 Torr between walls 2 cm apart, on 200 cells
across the gap, with c_p = 520.3 J/(kg K) (per kg in `thermoProperties`),
mu = 2.09e-5 Pa s and Pr = 2/3.

| Case | Test | Analytic solution |
|---|---|---|
| `advectionDiffusion` | steady flow (`rhoSimpleFoam` equations) and the transport of a species by the gas: metastables in a uniform flow U = 1 m/s from the inlet (N = 1e16 m^-3) to the outlet (N = 1e10 m^-3), D = 4.7e-3 m2/s, Peclet number 4.25 | N = N0 + (NL - N0)(exp(Pe x/L) - 1)/(exp(Pe) - 1); error 6e-5 |
| `channelStartup` | transient momentum (`rhoPimpleFoam` equations): gas at rest accelerated by a uniform force f = 0.418 N/m3 along the walls (`bodyForce`), periodic along the flow | away from the walls U = f t/rho, 1.955e-5 m/s at 0.1 microseconds; error 5e-6 |

Run `./Allrun` in a case directory: it makes the mesh, runs `somaFoam` and
compares the result with the analytic solution (`verify.py`), printing the
largest error relative to the scale of the solution and PASS or FAIL.
`./caseClean` removes the results and the mesh.

Known limitations of the gas flow, from cases that do not yet pass:

- The steady SIMPLE solution of a closed or periodic domain (no outlet that
  fixes the pressure), and of a flow without inertia, does not converge
  (plane Poiseuille flow, a uniformly heated gas at rest).
- The transient solution takes one flow step per plasma step and is unstable
  at acoustic Courant numbers far above 1.
- The gas energy is not yet verified: a gas heated uniformly in a closed box
  heats 14 % more slowly than at constant volume.
- With transient flow, the pressure must not be under-relaxed (factor 1 for
  `p` in `relaxationFactors`), since `p.relax()` of `rhoPimpleFoam` uses the
  factor of `p`.
