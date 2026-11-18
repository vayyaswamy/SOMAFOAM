"""Compare the electrode potential of the runs of Allrun with the RC circuit.

    python3 verify.py runs [--plot figure.pdf]

The gap is empty, so the electrode charge is Cg V with the vacuum capacitance
Cg = eps0 A/L of the 2 cm gap of area 4 cm^2. With the source V0 behind R and
a series capacitor C, the electrode potential obeys

    R Cg dV/dt + V (1 + Cg/C) = V0,

whose backward Euler solution, which the boundary condition solves exactly
when the gap is empty, is V^{n+1} = (V0 + a V^n)/(1 + Cg/C + a) with
a = R Cg/dt, and whose exact solution is V0 C/(C + Cg) (1 - exp(-t/tau)),
tau = R Cg C/(C + Cg).
"""
import math
import os
import sys

EPS0, A, L, R, V0 = 8.85418782e-12, 4e-4, 0.02, 2e5, 10.0
Cg = EPS0 * A / L
runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
cases = {"R_dt5e-10": (5e-10, 1e20), "R_dt1e-8": (1e-8, 1e20), "RC_dt5e-10": (5e-10, 1.770837564e-13)}

ok = True
curves = {}
errs = {}
for name, (dt, C) in cases.items():
    rows = []
    for line in open(os.path.join(runs, name, "electrodes", "0", "electrode.dat")):
        if line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    a = R * Cg / dt
    tau = R * Cg * C / (C + Cg)
    Vinf = V0 * C / (C + Cg)
    v, errBE, errExact = 0.0, 0.0, 0.0
    for k, r in enumerate(rows):
        if r[0] > 0:
            v = (V0 + a * v) / (1 + Cg / C + a)
        errBE = max(errBE, abs(r[1] - v) / V0)
        errExact = max(errExact, abs(r[1] - Vinf * (1 - math.exp(-r[0] / tau))) / V0)
    good = errBE < 1e-3 and len(rows) > 10
    ok &= good
    curves[name] = (rows, tau, Vinf)
    errs[name] = (dt, errExact)
    print(f"{name:11s} dt/tau = {dt/tau:6.3f}  steps {len(rows):4d}  V(end) = {rows[-1][1]:.6f} V"
          f"  max error against backward Euler {errBE:.2e}, against exp {errExact:.2e}"
          f"  {'ok' if good else 'too large'}")

# The error against the exact exponential falls with the step as backward
# Euler does, at first order
(d1, e1), (d2, e2) = errs["R_dt1e-8"], errs["R_dt5e-10"]
order = math.log(e1 / e2) / math.log(d1 / d2)
ok &= order > 0.9
print(f"order of the error against the exponential between dt = {d1:g} and {d2:g}: {order:.2f}")

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    for name, ls in (("R_dt5e-10", "o"), ("R_dt1e-8", "s"), ("RC_dt5e-10", "^")):
        rows, tau, Vinf = curves[name]
        t = [r[0] for r in rows]
        ax.plot([x * 1e9 for x in t], [Vinf * (1 - math.exp(-x / tau)) for x in t], "k-", lw=0.8)
        step = max(1, len(rows) // 20)
        ax.plot([r[0] * 1e9 for r in rows[::step]], [r[1] for r in rows[::step]], ls, ms=3, mfc="none", label=name.replace("_", ", "))
    ax.set_xlabel("t (ns)")
    ax.set_ylabel("electrode potential (V)")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
