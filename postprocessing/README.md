# Post-processing scripts

Python scripts that read SOMAFOAM cases directly (no ParaView or
reconstruction needed) and plot the usual quantities. They need Python 3 with
`numpy` and `matplotlib`, and cases written in ASCII (`writeFormat ascii;`).

Each script shows the figure, or saves it with `--output file.png`; `--help`
lists all options. Times are given in seconds or as `latest`, and the nearest
written time is used.

| Script | For |
|---|---|
| `plot_1d.py` | profiles and time histories of 1D cases |
| `plot_2d.py` | colour maps and line cuts of 2D cases, including adaptively refined meshes and parallel cases that have not been reconstructed |
| `plot_electrodes.py` | electrode voltage and current (any dimension) |
| `somafoam_io.py` | the readers used by the three scripts; can be imported in your own scripts |

## 1D cases

Profiles of several fields at several times, with the metastable density on a
logarithmic axis:

```
postprocessing/plot_1d.py examples/plasma/100TorrArgonPlasma \
    --fields electron Arp1 Arm Te Phi E --times 1e-6 5e-6 latest --log Arm
```

Time histories (maximum, gap average and mid-gap value at every written
time):

```
postprocessing/plot_1d.py CASE --fields electron Arm Te --history --log Arm
```

Field names are the file names in the time directories: species densities
(`electron`, `Arp1`, `Arm`), `Te`, `Phi`, `E`, and the averages written by the
`fieldAverage` function object (`TeMean`, `N_electronMean`, ...). The fields
written at a time are snapshots at that instant of the applied voltage; the
`...Mean` fields are averages over the interval since the previous output if
`resetOnOutput true;` is set in `system/controlDict`, which are averages over
whole periods when `writeInterval` is a multiple of the period.

`--csv file.csv` also writes the plotted data.

## 2D cases

Maps of several fields, drawn cell by cell, with the mesh and the refinement
level:

```
postprocessing/plot_2d.py examples/plasma/2DAdaptiveMeshArgon \
    --fields electron Te cellLevel --log electron --mesh
```

A cut along a line from (x0, y0) to (x1, y1), in metres:

```
postprocessing/plot_2d.py CASE --fields electron Te --line 0 4e-4 2e-3 4e-4
```

The mesh in use at the selected time is read from that time directory when
the mesh is refined during the run. For a vector field the magnitude is
plotted; `--component 0` selects a component.

## Electrode voltage and current

Requires the `electrodeVoltageCurrent` function object in
`system/controlDict` (see the main README).

```
postprocessing/plot_electrodes.py CASE --patch electrode --frequency 40e6
```

The top row shows the waveforms over the last periods of the run. The bottom
row has one point per period over the whole run: the RMS of the alternating
current, the period-averaged current and the period-averaged power. A
discharge that has reached a periodic steady state gives flat lines there.
`--csv` writes these per-period values.

## Using the readers directly

```python
import sys; sys.path.insert(0, "postprocessing")
import somafoam_io as sio

case = sio.Case("examples/plasma/100TorrArgonPlasma")
t = case.latest_time()
mesh = case.mesh(t)                 # points, faces, cell centres
ne = case.field("electron", t)      # numpy array, one value per cell
vi = sio.read_electrode(case.path, "electrode")   # t, V, I_total, ...
```
