# ArgonGlobal

Zero-dimensional argon model at 1 Torr on one cell: a 30 Td pulse of 2 us and
its afterglow. The drive and the options are explained in the comments of
`constant/plasmaProperties`.

```
blockMesh
somaFoam
```

The densities and Te of the cell are written every 10 steps to `cellProbe/0`.
