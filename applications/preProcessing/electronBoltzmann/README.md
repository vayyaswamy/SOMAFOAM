# electronBoltzmann: electron transport and rate tables from cross sections

`electronBoltzmann` solves the Boltzmann equation for electrons in a uniform
field (two-term approximation, or a Monte Carlo swarm) for a list of reduced
fields E/N. It writes the swarm data and the tables that `somaFoam` reads:
electron mobility, diffusion and rate coefficients as functions of the
electron temperature. The physics and every dictionary entry are described in
the header of `electronBoltzmann.C`.

## 1. Get the cross sections from LXCat

LXCat data may not be redistributed, so download them yourself:

1. Go to https://www.lxcat.net, open the data center and choose
   **Scattering cross sections**, electrons. The download page has a
   human-verification check, so use a browser.
2. Choose a database and the targets. Take a *complete set* (marked
   "complete set" by the database): it has the momentum transfer
   (ELASTIC or EFFECTIVE) and all inelastic processes, which the solver needs
   to get the energy distribution right. Examples:
   - Ar: Phelps database (3 processes), or TRINITI (4s levels resolved).
   - Air: Phelps N2 and O2, all processes marked "complete set".
3. Download as text and save it in the case, e.g.
   `constant/crossSections/Ar_Phelps.txt`.
4. Cite the database as LXCat asks (name, www.lxcat.net, retrieval date) and
   the original references listed in the file.

A file may hold several databases or targets. Only the targets listed in
`gases` are used, in the order of the file. Processes are numbered in that
order (the numbers are printed at the start of the run); they name the
`rate_<n>` tables and the `k<n>` columns of `transport.dat`.

## 2. Make a case

A case needs only `system/controlDict` and `system/boltzmannDict`; no mesh.
`examples/electronBoltzmann/Argon` is a ready case:

```
examples/electronBoltzmann/Argon
├── system/controlDict        time settings (not used)
├── system/boltzmannDict      cross-section file, gases, E/N range, method
└── constant/crossSections/   put the LXCat file here
```

The main entries of `system/boltzmannDict`:

| Entry | Meaning |
|---|---|
| `file` | LXCat file, relative to the case or absolute |
| `gases` | `((target fraction) ...)`, e.g. `((N2 0.79) (O2 0.21))` |
| `gasTemperature` | gas temperature [K] |
| `reducedField` | `{ min; max; n; }` logarithmic [Td], or a list `(1 10 100)` |
| `nEnergyCells` | energy grid (400–800); the upper end adapts to the distribution |
| `method` | `twoTerm` (fast) or `monteCarlo` (with `nElectrons`, `sampleTimes`, ...) |
| `nThreads` | reduced fields solved in parallel |
| `tables` | optional: write somaFoam tables under given names (see the header) |

## 3. Run

```
source <path-to>/SOMAFOAM/etc/bashrc
cd examples/electronBoltzmann/Argon
cp <your LXCat file> constant/crossSections/Ar_Phelps.txt
electronBoltzmann > log.electronBoltzmann
```

The log lists the processes that were read with their numbers, then one
line per reduced field. A two-term run over 121 fields takes about a minute.

## 4. Output (`constant/boltzmann/`)

| File | Content |
|---|---|
| `transport.dat` | per E/N [Td]: mean energy [eV], Te = 2/3 mean energy [K], mu N and D N, energy mobility and diffusion times N, and the rate coefficient k<n> [m3/s] of every process |
| `mu_electron`, `D_electron` | (Te [K], mu N or D N): mobility and diffusion tables of somaFoam |
| `rate_<n>` | (Te [K], ln k [m3/(kmol s)]) of process n: rate tables of somaFoam |
| `distribution_<E/N>.dat` | energy [eV] and F0 [eV^-3/2], normalised to int sqrt(eps) F0 = 1 |

Check the result before using it: plot the drift velocity mu N (E/N), the
mean energy and the Townsend coefficient alpha/N = k_ion/(mu N E/N) against
measured swarm data (LXCat also has swarm data).

## 5. Use the tables in somaFoam

Either

- copy `mu_electron`, `D_electron` and the `rate_<n>` files you need into the
  case `constant/` and name them in `constant/reactions`
  (`electronImpactInelasticTabularReaction ... ("$FOAM_CASE/constant/rate_3")`)
  and `constant/plasmaProperties` (`mobilityModel eTemp;` reads
  `constant/mu_electron`); or
- put `system/boltzmannDict` in the somaFoam case with a `tables { ... }`
  entry and set `electronBoltzmann yes;` in `constant/plasmaProperties`:
  somaFoam then runs the solver at start-up and writes the named tables
  before reading them.

A rate for a target that is not part of the gas (e.g. ionization of a
metastable, e + Ar* -> 2e + Ar+) is not part of the Boltzmann equation of the
background gas; compute it from the distributions:
k = sqrt(2e/m) * int eps sigma(eps) F0(eps) d eps.
