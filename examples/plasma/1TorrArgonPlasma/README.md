# 1 Torr argon: default schemes and Scharfetter-Gummel fluxes

Argon at 1 Torr in a 2 cm gap, 100 V at 40 MHz, one-dimensional. The same
discharge with each flux scheme (`fluxScheme` in `constant/plasmaProperties`):

| Case | `fluxScheme` | Mesh |
|---|---|---|
| `defaultSchemes` | `fvSchemes` (default) | 800 cells |
| `scharfetterGummel` | `scharfetterGummel` | 200 cells |

Run `somaFoam` in either directory (3 microseconds; about 20 and 55 minutes).
The default schemes need the finer mesh: with Scharfetter-Gummel fluxes 200
and 400 cells agree within 0.3 %, with the default schemes they differ by
15 %.
