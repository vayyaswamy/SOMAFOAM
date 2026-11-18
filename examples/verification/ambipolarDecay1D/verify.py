"""Ambipolar decay of the fundamental mode, against the analytic solution.

    python3 verify.py run [--plot figure.pdf]

Electrons (mu_e = 30 m2/(V s), Te = 1 eV) and ions (mu_i = 1.44 m2/(V s),
Ti = 300 K) with Einstein diffusion between walls 2 cm apart where both
densities vanish. In the quasi-neutral limit both decay as
n0 sin(pi x/L) exp(-t/tau_a), tau_a = L^2/(pi^2 D_a), with the ambipolar
diffusion coefficient D_a = (mu_i D_e + mu_e D_i)/(mu_i + mu_e). The Debye
length at n0 = 1e16 m^-3 is 7e-5 m, 0.4 % of the gap.
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

run = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None


def internal(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    if not m:
        return [float(re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)", s).group(1))] * NX
    return [float(v) for v in s[m.end():].split(")")[0].split()][:int(m.group(1))]


times = sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d)), key=float)
x = [(i + 0.5) * L / NX for i in range(NX)]
data = []
for t in times:
    ne = internal(os.path.join(run, t, "electron"))
    ni = internal(os.path.join(run, t, "Arp1"))
    data.append((float(t), sum(ne) * L / NX, sum(ni) * L / NX, ne, ni))

# Decay rate of the ion inventory over the second half of the run
half = [d for d in data if d[0] >= data[-1][0] / 2]
rate = math.log(half[0][2] / half[-1][2]) / (half[-1][0] - half[0][0])
errRate = abs(rate * tau - 1)

# Profile at the end against the sine of the same inventory
t, Ne, Ni, ne, ni = data[-1]
amp = Ni / (2 * L / math.pi)
errProfile = max(abs(n - amp * math.sin(math.pi * xi / L)) for n, xi in zip(ni, x)) / amp
quasi = max(abs(a - b) for a, b in zip(ne, ni)) / amp

print(f"D_a = {Da:.5f} m2/s, tau_a = {tau:.4e} s")
print(f"decay rate of the ions {rate:.5e} 1/s, exact {1/tau:.5e} 1/s, relative error {errRate:.2e}")
print(f"profile at t = {t:.3e} s: largest deviation from the sine {errProfile:.2e} of its peak, "
      f"largest |n_e - n_i| {quasi:.2e} of the peak")
ok = errRate < 1e-2 and errProfile < 2e-2

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    a1.semilogy([d[0] / tau for d in data], [d[2] / data[0][2] for d in data], "o", ms=3, mfc="none", label="ions")
    a1.semilogy([d[0] / tau for d in data], [d[1] / data[0][1] for d in data], "x", ms=3, label="electrons")
    a1.semilogy([d[0] / tau for d in data], [math.exp(-d[0] / tau) for d in data], "k-", lw=0.8, label="exp($-t/\\tau_a$)")
    a1.set_xlabel("$t/\\tau_a$")
    a1.set_ylabel("inventory / initial")
    a1.legend(frameon=False)
    for d in data[::4]:
        a2.plot([v * 1e3 for v in x[::5]], [n / data[0][4][NX // 2] for n in d[4][::5]], "o", ms=2.5, mfc="none")
        a2.plot([v * 1e3 for v in x], [math.exp(-d[0] / tau) * math.sin(math.pi * xi / L) for xi in x], "k-", lw=0.6)
    a2.set_xlabel("x (mm)")
    a2.set_ylabel("$n_\\mathrm{i}/n_0$")
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
