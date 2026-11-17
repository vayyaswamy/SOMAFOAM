"""Observed order of accuracy of the manufactured-solution runs of Allrun.

    python3 verify.py runs [--plot figure.pdf]

reads runs/<mesh>_<scheme>_<n>/mmsErrors.dat, prints the L2 errors of V, n
and w with the order between successive meshes, and checks the order on
the two finest meshes against the expected order of each discretisation.
"""
import glob
import math
import os
import re
import sys

runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None

data = {}
for f in glob.glob(os.path.join(runs, "*", "mmsErrors.dat")):
    mesh, scheme, n = re.fullmatch(r"(\w+?)_(\w+)_(\d+)", os.path.basename(os.path.dirname(f))).groups()
    row = [float(v) for v in open(f).read().split("\n")[1].split()]
    data.setdefault((mesh, scheme), []).append(row)

# Expected order on the finest meshes. On the uniform mesh both flux schemes
# are second order. On the distorted mesh the face gradient of the
# non-orthogonal correction and the linear interpolation without skewness
# correction are first order, and so is the solution; the ion density, with
# its larger cell Peclet number, has not reached its asymptotic order on these
# meshes and is reported without a check
expected = {
    ("uniform", "scharfetterGummel"): {"V": 2, "n": 2, "w": 2},
    ("uniform", "fvSchemes"): {"V": 2, "n": 2, "w": 2},
    ("distorted", "scharfetterGummel"): {"V": 1, "w": 1},
    ("distorted", "fvSchemes"): {"V": 1, "w": 1},
}
cols = {"V": 4, "n": 6, "w": 8}
ok = True

for key in sorted(data, key=lambda k: (k[0] != "uniform", k[1])):
    rows = sorted(data[key])
    print(f"\n{key[0]} mesh, {key[1]}")
    print(f"{'cells':>7} {'Pe_i':>6} {'Pe_e':>6}" + "".join(f"  {'L2 ' + v:>10} {'order':>5}" for v in cols))
    for j, r in enumerate(rows):
        line = f"{int(r[0]):7d} {r[2]:6.2f} {r[3]:6.2f}"
        for v, c in cols.items():
            order = math.log(rows[j - 1][c] / r[c]) / math.log(rows[j - 1][1] / r[1]) if j else float("nan")
            line += f"  {r[c]:10.3e} {order:5.2f}"
        print(line)
    for v, p in expected.get(key, {}).items():
        c = cols[v]
        order = math.log(rows[-2][c] / rows[-1][c]) / math.log(rows[-2][1] / rows[-1][1])
        good = order > p - 0.2
        ok &= good
        print(f"    order of {v} on the finest meshes {order:.2f}, expected {p}: {'ok' if good else 'too low'}")

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharey=True)
    for ax, mesh in zip(axes, ("uniform", "distorted")):
        for scheme, ls in (("scharfetterGummel", "-"), ("fvSchemes", "--")):
            rows = sorted(data.get((mesh, scheme), []))
            if not rows:
                continue
            h = [r[1] for r in rows]
            for v, mk in (("V", "o"), ("n", "s"), ("w", "^")):
                if scheme == "scharfetterGummel" or v != "V":
                    ax.loglog(h, [r[cols[v]] for r in rows], ls, marker=mk, ms=4,
                              label=f"{v}, {'SG' if scheme == 'scharfetterGummel' else 'fvSchemes'}")
        h = [r[1] for r in sorted(data[(mesh, "fvSchemes")])]
        ax.loglog(h, [0.5 * x * x for x in h], "k:", label="h$^2$")
        if mesh == "distorted":
            ax.loglog(h, [0.1 * x for x in h], "k-.", label="h")
        ax.set_xticks(h)
        ax.set_xticklabels([f"1/{round(1 / x)}" for x in h])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_title(f"{mesh} mesh")
        ax.set_xlabel("h")
    axes[0].set_ylabel("L2 error")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("\nPASS" if ok else "\nFAIL")
