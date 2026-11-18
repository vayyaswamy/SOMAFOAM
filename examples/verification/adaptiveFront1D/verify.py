"""Adaptive and uniform runs of the moving pulse against the exact solution.

    python3 verify.py runs [--plot figure.pdf]

The electron density of a Gaussian pulse of width s0 = 0.4 mm drifting at
v = -mu E, spreading with D and growing at nu, as in townsendAvalanche1D. For
each run, from the cell centres written at each time, the rates of growth,
drift and spreading of the moments, the largest error of the final profile
relative to its exact peak, and the number of cells.
"""
import math
import os
import re
import sys

KB_E, NA, KB = 8.617333262e-5, 6.02214076e26, 1.380649e-23
L = 0.01
mu, Te, Efield = 30.0, 1160.95, 1000.0
N = 133.3224 / (KB * 300.0)
nu = 1.870915e11 / NA * N
v = -mu * Efield
D = mu * KB_E * Te
x0, s0, n0 = 7e-3, 4e-4, 1e10


def internal(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    return [float(x) for x in s[m.end():].split(")")[0].split()][:int(m.group(1))]


def exact(x, t):
    s2 = s0 * s0 + 2 * D * t
    return n0 * s0 / math.sqrt(s2) * math.exp(nu * t) * math.exp(-(x - x0 - v * t)**2 / (2 * s2))


def slope(xs, ys):
    xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((a - xm) * (b - ym) for a, b in zip(xs, ys)) / sum((a - xm)**2 for a in xs)


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
res = {}
for name in ("uniform125", "uniform1000", "adaptive125"):
    run = os.path.join(runs, name)
    data = []
    for t in sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d)), key=float):
        xc = internal(os.path.join(run, t, "ccx"))
        n = [max(a - 1e4, 0.0) for a in internal(os.path.join(run, t, "electron"))]
        order = sorted(range(len(xc)), key=lambda i: xc[i])
        xc, n = [xc[i] for i in order], [n[i] for i in order]
        faces = [0.0] + [0.5 * (a + b) for a, b in zip(xc, xc[1:])] + [L]
        w = [b - a for a, b in zip(faces, faces[1:])]
        m0 = sum(a * b for a, b in zip(n, w))
        m1 = sum(a * b * c for a, b, c in zip(n, w, xc)) / m0
        var = sum(a * b * (c - m1)**2 for a, b, c in zip(n, w, xc)) / m0
        data.append((float(t), m0, m1, var, xc, n))
    t = [d[0] for d in data]
    eg = slope(t, [math.log(d[1]) for d in data]) / nu - 1
    ev = slope(t, [d[2] for d in data]) / v - 1
    ed = slope(t, [d[3] for d in data]) / 2 / D - 1
    tf, _, _, _, xc, n = data[-1]
    peak = exact(x0 + v * tf, tf)
    ep = max(abs(a - exact(b, tf)) for a, b in zip(n, xc)) / peak
    res[name] = (eg, ev, ed, ep, len(xc), data)
    print(f"{name:12s} {len(xc):5d} cells at the end: growth {eg:+.2e}, drift {ev:+.2e}, "
          f"diffusion {ed:+.2e}, profile {ep:.2e}")

a, u = res["adaptive125"], res["uniform1000"]
ok = a[3] < 2 * u[3] + 1e-3 and a[4] < u[4] and res["uniform125"][3] > a[3]
print(f"adaptive: profile error {a[3]:.2e} with {a[4]} cells against {u[3]:.2e} with {u[4]} uniform cells")

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    data = res["adaptive125"][5]
    for tt, m0, m1, var, xc, n in data:
        xs = [i * L / 1000 for i in range(1001)]
        a1.plot([x * 1e3 for x in xs], [exact(x, tt) / n0 for x in xs], "k-", lw=0.6)
        a1.plot([x * 1e3 for x in xc], [b / n0 for b in n], ".", ms=2)
    a1.set_xlabel("x (mm)")
    a1.set_ylabel("$n_\\mathrm{e}/n_0$")
    a1.set_title("(a) adaptive mesh, every 20 ns", fontsize=8)
    tt, m0, m1, var, xc, n = data[-1]
    widths = [b - a for a, b in zip([0.0] + [0.5 * (p + q) for p, q in zip(xc, xc[1:])], [0.5 * (p + q) for p, q in zip(xc, xc[1:])] + [L])]
    a2.semilogy([x * 1e3 for x in xc], [w_ * 1e6 for w_ in widths], "-")
    a2.set_xlabel("x (mm)")
    a2.set_ylabel("cell size ($\\mu$m)")
    a2.set_title(f"(b) mesh at 100 ns, {len(xc)} cells", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
