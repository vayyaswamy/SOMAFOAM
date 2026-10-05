# 1 Torr argon: default schemes and Scharfetter-Gummel fluxes

The same discharge solved with the two flux schemes that can be selected with
`fluxScheme` in `constant/plasmaProperties`:

| Case | `fluxScheme` | Mesh |
|---|---|---|
| `defaultSchemes` | `fvSchemes` (the default) | 800 uniform cells (25 um) |
| `scharfetterGummel` | `scharfetterGummel` | 200 uniform cells (100 um) |

Argon at 1 Torr and 300 K between a driven electrode (x = 0,
100 V sin(2 pi 40 MHz t), no DC bias) and a grounded electrode (x = 2 cm),
one-dimensional. Electrons, Ar+ and the metastable Arm with the `mixed` model
(drift-diffusion), electron temperature from `efullImplicit`, secondary
emission coefficient 0.05, fixed time step 1e-10 s, 3 microseconds (30000
steps) with output every 0.25 microseconds. The two cases differ only in
`fluxScheme` and in the number of cells.

Run (serial; the mesh is included, `blockMesh` recreates it):

```
cd scharfetterGummel      # or defaultSchemes
somaFoam
```

About 20 minutes of CPU time for `scharfetterGummel` and 55 minutes for
`defaultSchemes`. The discharge is still developing at 3 microseconds; raise
`endTime` for longer runs.

## Why the meshes differ

The sheaths are about 2 mm wide and the electron density falls by four orders
of magnitude across them. With the default schemes the result depends
strongly on the cell size; with Scharfetter-Gummel fluxes 200 cells are
enough. At 3 microseconds:

| Scheme | Cells | Peak electron density (1/m3) | RF current, RMS (mA) | Power (W) | Peak Te in the sheath (K) |
|---|---|---|---|---|---|
| default | 200 | 2.91e16 | 0.255 | 0.00229 | 79100 |
| default | 400 | 3.42e16 | 0.297 | 0.00259 | 64200 |
| default | 800 | 3.88e16 | 0.336 | 0.00292 | 55300 |
| Scharfetter-Gummel | 200 | 4.21e16 | 0.357 | 0.00316 | 51000 |
| Scharfetter-Gummel | 400 | 4.21e16 | 0.358 | 0.00315 | 52300 |

The Scharfetter-Gummel results on 200 and 400 cells agree within 0.3 %. The
default results are still rising with the number of cells at 800 cells and
are about 7 % lower in density there; which value they converge to has not
been established.

## Plots

```
python3 ../../../../postprocessing/plot_1d.py . --fields electron Arp1 Te Phi --log electron Arp1
python3 ../../../../postprocessing/plot_electrodes.py . --frequency 40e6
```
