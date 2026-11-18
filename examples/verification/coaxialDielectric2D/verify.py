"""Potential of the coaxialDielectric2D runs against the exact solution.

    python3 verify.py runs [--plot figure.pdf]

With the electrode at r = a and potential V0, a dielectric of permittivity
eps_d to r = b and gas of permittivity eps to the grounded wall at r = c, and
no space charge, the charge per length is
q = 2 pi V0/(ln(b/a)/eps_d + ln(c/b)/eps), and the potential is
V0 - q/(2 pi eps_d) ln(r/a) in the dielectric and q/(2 pi eps) ln(c/r) in
the gas. The error is the largest |Phi - Phi_exact| over the cells of both
regions relative to V0.
"""
import glob
import math
import os
import re
import sys

a, b, c, V0 = 1e-3, 2e-3, 5e-3, 10.0
eps, epsd = 8.85e-12, 4.425e-11
q = 2 * math.pi * V0 / (math.log(b / a) / epsd + math.log(c / b) / eps)


def exact(r):
    if r <= b:
        return V0 - q / (2 * math.pi * epsd) * math.log(r / a)
    return q / (2 * math.pi * eps) * math.log(c / r)


def field(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    body = s[m.end():]
    if m.group(1) == "vector":
        return [tuple(float(x) for x in v.split()) for v in re.findall(r"\(([^()]*)\)", body)[:int(m.group(2))]]
    return [float(x) for x in body.split(")")[0].split()][:int(m.group(2))]


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
results = []
for run in sorted(glob.glob(os.path.join(runs, "n_*")), key=lambda d: int(d.split("_")[-1])):
    nt = int(run.split("_")[-1])
    nd, ng = nt // 2, 3 * nt // 2
    t = max((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d)), key=float)
    errs, prof = [], []
    for region, r1, r2, nr in (("dielectric_1", a, b, nd), ("", b, c, ng)):
        phi = field(os.path.join(run, t, region, "Phi"))
        # cell centres written by writeCellCentres at time 0
        cx = field(os.path.join(run, "0", region, "ccx"))
        cy = field(os.path.join(run, "0", region, "ccy"))
        for k, val in enumerate(phi):
            rc = math.hypot(cx[k], cy[k])
            errs.append(abs(val - exact(rc)) / V0)
            if k < nr:
                prof.append((rc, val))
    err = max(errs)
    results.append((nt, err, sorted(prof)))
    print(f"{nt:3d} cells along the arc: max |Phi - Phi_exact|/V0 = {err:.3e}")

ok = len(results) == 3
for (n1, e1, _), (n2, e2, _) in zip(results, results[1:]):
    p = math.log(e1 / e2) / math.log(n2 / n1)
    print(f"order between {n1} and {n2}: {p:.2f}")
    ok &= p > 1.7
ok &= results[-1][1] < 1e-3

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    rr = [a + i * (c - a) / 400 for i in range(401)]
    a1.plot([r * 1e3 for r in rr], [exact(r) for r in rr], "k-", lw=0.8, label="exact")
    prof = results[-1][2]
    a1.plot([r * 1e3 for r, _ in prof[::3]], [v for _, v in prof[::3]], "o", ms=3, mfc="none", label="SOMAFOAM")
    a1.axvspan(a * 1e3, b * 1e3, color="0.9")
    a1.set_xlabel("r (mm)")
    a1.set_ylabel("potential (V)")
    a1.legend(frameon=False)
    a1.set_title("(a) potential, dielectric shaded", fontsize=8)
    a2.loglog([r[0] for r in results], [r[1] for r in results], "o-", ms=4, label="max error")
    a2.loglog([r[0] for r in results], [results[0][1] * (results[0][0] / r[0])**2 for r in results], "k:", label="$N^{-2}$")
    a2.set_xlabel("cells along the arc")
    a2.set_ylabel("error / $V_0$")
    a2.legend(frameon=False)
    a2.set_title("(b) convergence", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
