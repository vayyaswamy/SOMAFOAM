"""Electron avalanche in a uniform field against the exact solution.

    python3 verify.py runs [--plot figure.pdf]

With the mobility mu = 30 m2/(V s), Te = 0.1 eV (D = mu kB Te/e), the field
E = 10 V/1 cm and the ionization frequency nu = k N of the single reaction,
the electron density obeys dn/dt + d(-mu E n)/dx - D d2n/dx2 = nu n, whose
solution for a Gaussian start is a Gaussian whose content grows as
exp(nu t), whose centre moves at -mu E and whose variance grows as 2 D t.
The three are measured from the moments of the computed density above the
floor of 1e4 m^-3 at the write times and compared by their rates. The
implicit Euler step adds the numerical diffusion v^2 dt/2 to D, so the
spreading is checked at two time steps for an excess that falls with the step.
"""
import math
import os
import re
import sys

KB_E, NA, KB = 8.617333262e-5, 6.02214076e26, 1.380649e-23
L, NX = 0.01, 1000
mu, Te, Efield = 30.0, 1160.95, 1000.0
N = 133.3224 / (KB * 300.0)
nu = 1.870915e11 / NA * N
v = -mu * Efield
D = mu * KB_E * Te


def internal(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    return [float(x) for x in s[m.end():].split(")")[0].split()][:int(m.group(1))]


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
x = [(i + 0.5) * L / NX for i in range(NX)]


def slope(xs, ys):
    xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((a - xm) * (b - ym) for a, b in zip(xs, ys)) / sum((a - xm)**2 for a in xs)


results = []
for run in sorted((os.path.join(runs, d) for d in os.listdir(runs) if d.startswith("dt_")), key=lambda d: -float(d.split("_")[-1])):
    dt = float(run.split("_")[-1])
    data = []
    for t in sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d)), key=float):
        n = [max(a - 1e4, 0.0) for a in internal(os.path.join(run, t, "electron"))]
        m0 = sum(n)
        m1 = sum(a * b for a, b in zip(n, x)) / m0
        var = sum(a * (b - m1)**2 for a, b in zip(n, x)) / m0
        data.append((float(t), m0, m1, var, n))
    t = [d[0] for d in data]
    growth = slope(t, [math.log(d[1]) for d in data])
    drift = slope(t, [d[2] for d in data])
    spread = slope(t, [d[3] for d in data]) / 2
    eg, ev, ed = growth / nu - 1, drift / v - 1, spread / D - 1
    results.append((dt, eg, ev, ed, data))
    print(f"dt = {dt:g} s: ionization frequency {growth:.6e} 1/s ({eg:+.2e}), drift velocity {drift:.6e} m/s ({ev:+.2e}), "
          f"diffusion coefficient {spread:.5e} m2/s ({ed:+.2e}; numerical diffusion v^2 dt/2 is {v*v*dt/2/D:.2e} of D)")
print(f"exact: ionization frequency {nu:.6e} 1/s, drift velocity {v:.6e} m/s, diffusion coefficient {D:.6e} m2/s")

(d1, g1, v1, e1, _), (d2, g2, v2, e2, data) = results
order = math.log(e1 / e2) / math.log(d1 / d2)
print(f"order of the excess diffusion in the time step: {order:.2f}")
ok = max(abs(g2), abs(v2)) < 2e-3 and abs(e2) < 6e-3 and order > 0.8

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    x0, s0 = 7e-3, 2e-4
    for tt, m0, m1, var, n in data:
        s2 = s0 * s0 + 2 * D * tt
        ex = [1e10 * s0 / math.sqrt(s2) * math.exp(nu * tt) * math.exp(-(xi - x0 - v * tt)**2 / (2 * s2)) for xi in x]
        ax.plot([xi * 1e3 for xi in x], [e / 1e10 for e in ex], "k-", lw=0.7)
        ax.plot([xi * 1e3 for xi in x[::25]], [a / 1e10 for a in n[::25]], "o", ms=2.5, mfc="none")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("$n_\\mathrm{e}/n_0$")
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
