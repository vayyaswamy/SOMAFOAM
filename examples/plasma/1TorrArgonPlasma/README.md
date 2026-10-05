# 1 Torr argon: default schemes and Scharfetter-Gummel fluxes

Argon at 1 Torr in a 2 cm gap, 100 V at 40 MHz, one-dimensional, 200 cells.
The same discharge with each flux scheme (`fluxScheme` in
`constant/plasmaProperties`); the two cases differ only in that entry:

| Case | `fluxScheme` |
|---|---|
| `defaultSchemes` | `fvSchemes` (default) |
| `scharfetterGummel` | `scharfetterGummel` |

Run `somaFoam` in either directory (3 microseconds, about 20 minutes). The
two agree within about 1 %.
