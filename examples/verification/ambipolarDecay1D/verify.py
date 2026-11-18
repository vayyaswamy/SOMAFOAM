"""Ambipolar decay of the fundamental mode, against the analytic solution.

    python3 verify.py runs [--plot figure.pdf]

Electrons (mu_e = 30 m2/(V s), Te = 1 eV) and ions (mu_i = 1.44 m2/(V s),
Ti = 300 K) with Einstein diffusion between absorbing walls 2 cm apart,
with the wall conditions of the plasma examples and no secondary emission.
In the quasi-neutral limit, where the densities vanish at the walls, both
decay as
n0 sin(pi x/L) exp(-t/tau_a), tau_a = L^2/(pi^2 D_a), with the ambipolar
diffusion coefficient D_a = (mu_i D_e + mu_e D_i)/(mu_i + mu_e). The Debye
length is 2.4e-5 m at n0 = 1e17 m^-3 and 2.4e-6 m at 1e19 m^-3, and the
sheaths of a few Debye lengths shorten the quasi-neutral region, which raises
the decay rate by an amount that falls with the Debye length. The runs at
three densities and, at the highest, on three meshes check that the error
falls with the Debye length and with the cell size, the latter because the
sheath at the highest density sits in the first cell, and is below 1.5 % on
the finest mesh.
"""
import math
import os
import re
import sys

KB_E = 8.617333262e-5          # kB/e [V/K]
L, NX = 0.02, 200
mue, mui, Te, Ti = 30.0, 1.44, 11605.0, 300.0
De, Di = mue * KB_E * Te, mui * KB_E * Ti
Da = (mui * De + mue * Di) / (mui + mue)
tau = L * L / (math.pi**2 * Da)

runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None


def internal(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    if not m:
        return [float(re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)", s).group(1))] * NX
    return [float(v) for v in s[m.end():].split(")")[0].split()][:int(m.group(1))]


x = [(i + 0.5) * L / NX for i in range(NX)]
print(f"D_a = {Da:.5f} m2/s, tau_a = {tau:.4e} s")
results = []
for name in sorted((d for d in os.listdir(runs) if re.fullmatch(r"n0_[0-9.e+]+_nx_[0-9]+", d)), key=lambda d: (int(d.split("_")[3]), float(d.split("_")[1]))):
    run = os.path.join(runs, name)
    n0, nx = float(name.split("_")[1]), int(name.split("_")[3])
    x = [(i + 0.5) * L / nx for i in range(nx)]
    times = sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d)), key=float)
    data = []
    for t in times:
        ne = internal(os.path.join(run, t, "electron"))
        ni = internal(os.path.join(run, t, "Arp1"))
        data.append((float(t), sum(ne) * L / nx, sum(ni) * L / nx, ne, ni))
    half = [d for d in data if d[0] >= data[-1][0] / 2]
    rate = math.log(half[0][2] / half[-1][2]) / (half[-1][0] - half[0][0])
    errRate = rate * tau - 1
    t, Ne, Ni, ne, ni = data[-1]
    amp = Ni / (2 * L / math.pi)
    errProfile = max(abs(n - amp * math.sin(math.pi * xi / L)) for n, xi in zip(ni, x)) / amp
    # Debye length at the mean density of the second half
    nmean = 0.5 * (max(half[0][4]) + max(half[-1][4]))
    lamD = math.sqrt(8.8541878e-12 * KB_E * Te / (nmean * 1.602176634e-19))
    results.append((n0, lamD, errRate, errProfile, data, nx))
    print(f"n0 = {n0:.0e}, {nx:4d} cells: Debye length {lamD:.2e} m ({lamD/L:.1e} of L), decay rate {rate:.5e} 1/s, "
          f"relative error {errRate:+.3e}, profile error {errProfile:.2e}")

dens = [r for r in results if r[5] == 200]
mesh = [r for r in results if r[0] == max(q[0] for q in results)]
ed = [abs(r[2]) for r in dens]
em = [abs(r[2]) for r in mesh]
ok = all(a > b for a, b in zip(ed, ed[1:])) and all(a > b for a, b in zip(em, em[1:])) and em[-1] < 1.5e-2
if len(mesh) > 1:
    p = math.log(em[0] / em[-1]) / math.log(mesh[-1][5] / mesh[0][5])
    print(f"order of the decay-rate error in the cell size at n0 = {mesh[0][0]:.0e}: {p:.2f}")
if len(mesh) > 2:
    d1, d2 = em[-3] - em[-2], em[-2] - em[-1]
    print(f"Richardson extrapolation of the error to a fine mesh: {em[-1] - d2*d2/(d1 - d2):.2e}")

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    data = mesh[-1][4]
    a1.semilogy([d[0] / tau for d in data], [d[2] / data[0][2] for d in data], "o", ms=3, mfc="none", label="ions")
    a1.semilogy([d[0] / tau for d in data], [d[1] / data[0][1] for d in data], "x", ms=3, label="electrons")
    a1.semilogy([d[0] / tau for d in data], [math.exp(-d[0] / tau) for d in data], "k-", lw=0.8, label="exp($-t/\\tau_a$)")
    a1.set_xlabel("$t/\\tau_a$")
    a1.set_ylabel("inventory / initial")
    a1.set_title(f"(a) $n_0$ = {mesh[-1][0]:.0e} m$^{{-3}}$, {mesh[-1][5]} cells", fontsize=9)
    a1.legend(frameon=False)
    a2.loglog([max(r[1], L / r[5]) / L for r in dens], ed, "o-", ms=4, label="200 cells, $n_0$ = 1e17 to 1e19")
    a2.loglog([L / r[5] / L for r in mesh], em, "s--", ms=4, label="$n_0$ = 1e19, 200 to 800 cells")
    a2.set_xlabel("max($\\lambda_\\mathrm{D}$, $\\Delta x$)/L")
    a2.set_ylabel("relative error of the decay rate")
    a2.set_title("(b) error of the decay rate", fontsize=9)
    a2.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
