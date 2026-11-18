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
three densities check that the error falls and is below 3 % at the highest.
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
for name in sorted(os.listdir(runs), key=lambda d: float(d.split("_")[1])):
    run = os.path.join(runs, name)
    n0 = float(name.split("_")[1])
    times = sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d)), key=float)
    data = []
    for t in times:
        ne = internal(os.path.join(run, t, "electron"))
        ni = internal(os.path.join(run, t, "Arp1"))
        data.append((float(t), sum(ne) * L / NX, sum(ni) * L / NX, ne, ni))
    half = [d for d in data if d[0] >= data[-1][0] / 2]
    rate = math.log(half[0][2] / half[-1][2]) / (half[-1][0] - half[0][0])
    errRate = rate * tau - 1
    t, Ne, Ni, ne, ni = data[-1]
    amp = Ni / (2 * L / math.pi)
    errProfile = max(abs(n - amp * math.sin(math.pi * xi / L)) for n, xi in zip(ni, x)) / amp
    # Debye length at the mean density of the second half
    nmean = 0.5 * (max(half[0][4]) + max(half[-1][4]))
    lamD = math.sqrt(8.8541878e-12 * KB_E * Te / (nmean * 1.602176634e-19))
    results.append((n0, lamD, errRate, errProfile, data))
    print(f"n0 = {n0:.0e}: Debye length {lamD:.2e} m ({lamD/L:.1e} of L), decay rate {rate:.5e} 1/s, "
          f"relative error {errRate:+.3e}, profile error {errProfile:.2e}")

errs = [abs(r[2]) for r in results]
ok = errs[-1] < 3e-2 and all(a > b for a, b in zip(errs, errs[1:]))
if len(results) > 1:
    p = math.log(errs[0] / errs[-1]) / math.log(results[0][1] / results[-1][1])
    print(f"order of the decay-rate error in the Debye length: {p:.2f}")

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    data = results[-1][4]
    a1.semilogy([d[0] / tau for d in data], [d[2] / data[0][2] for d in data], "o", ms=3, mfc="none", label="ions")
    a1.semilogy([d[0] / tau for d in data], [d[1] / data[0][1] for d in data], "x", ms=3, label="electrons")
    a1.semilogy([d[0] / tau for d in data], [math.exp(-d[0] / tau) for d in data], "k-", lw=0.8, label="exp($-t/\\tau_a$)")
    a1.set_xlabel("$t/\\tau_a$")
    a1.set_ylabel("inventory / initial")
    a1.set_title(f"(a) $n_0$ = {results[-1][0]:.0e} m$^{{-3}}$", fontsize=9)
    a1.legend(frameon=False)
    a2.loglog([r[1] / L for r in results], [abs(r[2]) for r in results], "o-", ms=4)
    a2.set_xlabel("$\\lambda_\\mathrm{D}/L$")
    a2.set_ylabel("relative error of the decay rate")
    a2.set_title("(b) error against the Debye length", fontsize=9)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
