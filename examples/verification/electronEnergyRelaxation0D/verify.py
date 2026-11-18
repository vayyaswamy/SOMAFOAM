"""Electron temperature of a zero-dimensional argon gas against the exact solution.

    python3 verify.py runs [--plot figure.pdf]

With a constant electron mobility mu, the collision frequency nu = e/(m mu)
of the muBased elastic loss, no reactions and a reduced field E/N, the
electron energy equation of the global model is

    d/dt (3/2 kB Te) = e mu E^2 - 3 kB (m/M) nu (Te - T),

so Te relaxes to T + mu^2 E^2 M/(3 kB) at the rate 2 e/(M mu) while the
field is on, and back to T after. The error is the largest |Te - Te_exact|
relative to the rise. The energy budget written at the write times
(writeEnergyBudget) must close: storage = heating + losses, with the
residual relative to the heating.
"""
import glob
import math
import os
import re
import sys

KB, E, U = 1.380649e-23, 1.602176634e-19, 1.66053906660e-27
p, T, mu, M = 133.3224, 300.0, 30.0, 39.948 * U
N = p / (KB * T)
EN, tOff = 2.75e-21, 2e-5
Efield = EN * N
lam = 2 * E / (M * mu)
dTinf = mu**2 * Efield**2 * M / (3 * KB)


def exact(t):
    if t <= tOff:
        return T + dTinf * (1 - math.exp(-lam * t))
    return T + dTinf * (1 - math.exp(-lam * tOff)) * math.exp(-lam * (t - tOff))


def probe(run):
    rows = []
    for f in glob.glob(os.path.join(run, "cellProbe", "*", "Te")):
        for line in open(f):
            if line.strip() and not line.startswith("#"):
                v = line.split()
                rows.append((float(v[0]), float(v[1])))
    return sorted(rows)


def uniform(path):
    s = open(path).read()
    m = re.search(r"internalField\s+uniform\s+(\S+);", s)
    if m:
        return float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\(", s)
    return float(s[m.end():].split(")")[0].split()[0])


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
print(f"N = {N:.4e} m^-3, E = {Efield:.4f} V/m, rise {dTinf:.2f} K, time constant {1/lam:.4e} s")

results = []
for run in sorted(glob.glob(os.path.join(runs, "dt_*")), key=lambda d: float(d.split("_")[-1])):
    dt = float(run.split("_")[-1])
    rows = probe(run)
    # the first step after the field is switched off is excluded from the
    # comparison: the switch falls inside it
    err = max(abs(te - exact(t)) for t, te in rows if abs(t - tOff) > 1.5 * dt) / dTinf
    peak = max(te for t, te in rows)
    # energy budget at the write times
    res, heat = [], []
    for d in os.listdir(run):
        if re.fullmatch(r"[0-9.eE+-]+", d) and float(d) > 0 and os.path.exists(os.path.join(run, d, "eResidual")):
            res.append(abs(uniform(os.path.join(run, d, "eResidual"))))
            heat.append(abs(uniform(os.path.join(run, d, "eHeating"))))
    budget = max(res) / max(heat) if heat else float("nan")
    results.append((dt, err, budget, rows))
    print(f"dt = {dt:g} s: {len(rows)} steps, peak Te {peak:.2f} K (exact {exact(tOff):.2f} K), "
          f"max error {err:.2e} of the rise, energy budget residual {budget:.1e} of the heating")

ok = len(results) == 2
(d1, e1, b1, _), (d2, e2, b2, _) = results
order = math.log(e2 / e1) / math.log(d2 / d1)
print(f"order of the error in the time step: {order:.2f}")
ok &= e1 < 2e-3 and order > 0.9 and max(b1, b2) < 1e-6

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    rows = results[0][3]
    tt = [i * 5e-5 / 500 for i in range(501)]
    ax.plot([t * 1e6 for t in tt], [exact(t) / 11604.5 for t in tt], "k-", lw=0.8, label="exact")
    step = max(1, len(rows) // 40)
    ax.plot([t * 1e6 for t, _ in rows[::step]], [te / 11604.5 for _, te in rows[::step]], "o", ms=3, mfc="none", label="SOMAFOAM")
    ax.set_xlabel("t ($\\mu$s)")
    ax.set_ylabel("$k_\\mathrm{B}T_\\mathrm{e}/e$ (V)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
