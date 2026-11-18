"""Compare the turnerCCP runs with the particle-in-cell benchmark.

    python3 compare.py runs [--plot figure.pdf]

From the field averages over the last 32 periods (the last write), the
mid-plane ion density and electron temperature and the ion current at the
electrodes, e times the averaged ion flux in the cells next to the two
electrodes, against table III of Turner et al. (2013), and for case 1 the ion
density profile against the benchmark profile (reference/). The fluid model
and the particle codes differ in their physics at low pressure, where the
electrons are non-local, so the differences are reported and not tested.
"""
import glob
import json
import math
import os
import re
import sys

E = 1.602176634e-19
KB_E = 8.617333262e-5
L = 0.067
here = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(os.path.join(here, "reference", "table3.json")))
prof1 = json.load(open(os.path.join(here, "reference", "case1_ion_density.json")))["ion_density_m3"]


def internal(path, vector=False):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", s)
    body = s[m.end():]
    if vector:
        return [float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", body)[:int(m.group(1))]]
    return [float(x) for x in body.split(")")[0].split()][:int(m.group(1))]


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
rows = []
print(f"{'case':>4} {'periods':>7} {'n_i (1e15 m^-3)':>22} {'Te (eV)':>18} {'J_i (A/m2)':>18}")
for run in sorted(glob.glob(os.path.join(runs, "case*"))):
    c = int(run[-1])
    times = sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d) and float(d) > 0), key=float)
    if not times or not os.path.exists(os.path.join(run, times[-1], "N_Hep1Mean")):
        print(f"{c:4d}  no averaged fields")
        continue
    t = times[-1]
    ni = internal(os.path.join(run, t, "N_Hep1Mean"))
    te = internal(os.path.join(run, t, "TeMean"))
    ji = internal(os.path.join(run, t, "J_Hep1Mean"), vector=True)
    n = len(ni)
    mid = lambda f: 0.5 * (f[n // 2 - 1] + f[n // 2])
    ni_m, te_m = mid(ni) / 1e15, mid(te) * KB_E
    Ji = 0.5 * E * (abs(ji[0]) + abs(ji[-1]))
    periods = round(float(t) * 13.56e6)
    r = (c, ni_m, te_m, Ji, ni)
    rows.append(r)
    k = c - 1
    print(f"{c:4d} {periods:7d} {ni_m:9.3f} ({ref['ni_1e15_m3'][k]:5.3f}, {ni_m/ref['ni_1e15_m3'][k]-1:+5.0%})"
          f" {te_m:6.2f} ({ref['Te_eV'][k]:4.2f}, {te_m/ref['Te_eV'][k]-1:+4.0%})"
          f" {Ji:6.3f} ({ref['Ji_A_m2'][k]:5.3f}, {Ji/ref['Ji_A_m2'][k]-1:+4.0%})")
    if c == 1:
        xr = [i * L / (len(prof1) - 1) for i in range(len(prof1))]
        xc = [(i + 0.5) * L / n for i in range(n)]
        interp = []
        for x in xr:
            j = min(max(int(x / (L / n) - 0.5), 0), n - 2)
            w = (x - xc[j]) / (xc[j + 1] - xc[j])
            interp.append((1 - w) * ni[j] + w * ni[j + 1])
        rms = math.sqrt(sum((a - b)**2 for a, b in zip(interp, prof1)) / len(prof1)) / max(prof1)
        print(f"     case 1 ion density profile: RMS difference {rms:.1%} of the benchmark peak")

if plot and rows:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    for c, ni_m, te_m, Ji, ni in rows:
        if c == 1:
            n = len(ni)
            a1.plot([(i + 0.5) * L / n * 100 for i in range(n)], [v / 1e14 for v in ni], "-", label="SOMAFOAM")
            a1.plot([i * L / (len(prof1) - 1) * 100 for i in range(len(prof1))][::4], [v / 1e14 for v in prof1][::4],
                    "ko", ms=3, mfc="none", label="PIC benchmark")
    a1.set_xlabel("x (cm)")
    a1.set_ylabel("$n_\\mathrm{i}$ ($10^{14}$ m$^{-3}$)")
    a1.set_title("(a) case 1, time-averaged ion density", fontsize=8)
    a1.legend(frameon=False, fontsize=7)
    cs = [r[0] for r in rows]
    a2.semilogy(cs, [r[1] for r in rows], "o-", label="$n_\\mathrm{i}$, SOMAFOAM")
    a2.semilogy([1, 2, 3, 4], ref["ni_1e15_m3"], "ks", mfc="none", label="$n_\\mathrm{i}$, PIC")
    a2.set_xticks([1, 2, 3, 4])
    a2.set_xlabel("benchmark case")
    a2.set_ylabel("mid-plane $n_\\mathrm{i}$ ($10^{15}$ m$^{-3}$)")
    a2.legend(frameon=False, fontsize=7)
    a2.set_title("(b) mid-plane ion density", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)
