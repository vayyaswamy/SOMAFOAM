# CF4O2

Zero-dimensional CF4/O2 model at 1 Torr and 300 K, 80 % CF4 and 20 % O2, on
one cell: 29 species and 125 reactions, driven by 150 Td for 3 us and then an
afterglow at 1 Td. The electron-impact rates of CF4 and O2 come from the LXCat
Hayashi and Phelps cross sections over the Boltzmann distributions, and those
of the radicals from Kimura and Noto (2006). `MECHANISM.md` lists every
reaction with its rate, energy loss and source; the choices are described in
the header of `mechanism/build_cf4o2_mechanism.py`.

```
blockMesh
somaFoam
```

The densities and Te of the cell are written every step to `cellProbe/0`.

To rebuild the chemistry (for another pressure or gas temperature, or other
cross sections), download the LXCat Hayashi CF4 and Phelps O2 sets into one
file, `mechanism/boltzmann/constant/crossSections/CF4_Hayashi_O2_Phelps.txt`
(LXCat data may not be redistributed), and run

```
cd mechanism/boltzmann && electronBoltzmann > log.electronBoltzmann && cd ../..
python3 mechanism/build_cf4o2_mechanism.py mechanism/boltzmann/constant/boltzmann . 133.322 300
```
