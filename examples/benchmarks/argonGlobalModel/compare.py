"""Steady state of the argon global model against Lieberman and Lichtenberg.

    python3 compare.py runs [--plot figure.pdf]

The global model of Lieberman and Lichtenberg (Principles of Plasma
Discharges and Materials Processing, 2nd ed., 2005, chapter 10) for an
electropositive argon discharge in a cylinder of radius R and length L:

    particle balance   N K_iz(Te) = u_B A_eff/V,  u_B = (kB Te/M)^1/2,
                       A_eff = 2 pi (R^2 h_L + R L h_R),
                       h_L = 0.86 (3 + L/(2 lambda_i))^-1/2,
                       h_R = 0.80 (4 + R/lambda_i)^-1/2,  lambda_i = 1/(N sigma_i)
    power balance      P/V = n_e [N (K_iz e_iz + K_ex e_ex + 3 (m/M) kB (Te - T) K_el)
                                  + (u_B A_eff/V) 2 kB Te],

i.e. the electron power with the collisional losses and the kinetic energy
2 kB Te of each electron lost to the walls. The rate coefficients are the Maxwellian fits of
Lieberman and Lichtenberg (makeRates.py), evaluated directly here and as
tables of 600 points in the case. The
particle balance gives Te alone, and the power balance then gives n_e.
The steady state of the runs is the last probe value.
"""
import glob
import math
import os
import re
import sys

KB, E, U, AV = 1.380649e-23, 1.602176634e-19, 1.66053906660e-27, 6.02214076e26
ME = 9.1093837e-31
M = 39.948 * U
R, Lc, sigma, Tg = 0.15, 0.3, 1e-18, 300.0
V = math.pi * R * R * Lc
EIZ, EEX = 2.5248e-18, 1.851912e-18
here = os.path.dirname(os.path.abspath(__file__))


sys.path.insert(0, here)
from makeRates import Kiz, Kex, Kel


def K(kind, Te):
    """Rate coefficient [m3/s] of the Maxwellian fits, Te in K."""
    return {"iz": Kiz, "ex": Kex, "el": Kel}[kind](Te / 11604.518)


def global_model(P, N):
    lam = 1 / (N * sigma)
    hL = 0.86 / math.sqrt(3 + Lc / (2 * lam))
    hR = 0.80 / math.sqrt(4 + R / lam)
    Aeff = 2 * math.pi * (R * R * hL + R * Lc * hR)
    kw = lambda Te: math.sqrt(KB * Te / M) * Aeff / V
    f = lambda Te: N * K("iz", Te) - kw(Te)
    a, b = 5000.0, 2e5
    for _ in range(200):
        m = 0.5 * (a + b)
        if f(a) * f(m) <= 0:
            b = m
        else:
            a = m
    Te = 0.5 * (a + b)
    loss = N * (K("iz", Te) * EIZ + K("ex", Te) * EEX + 3 * ME / M * KB * (Te - Tg) * K("el", Te)) + kw(Te) * 2 * KB * Te
    return Te, P / V / loss


def probe(run, field):
    f = glob.glob(os.path.join(run, "cellProbe", "*", field))[0]
    rows = [l.split() for l in open(f) if l.strip() and not l.startswith("#")]
    return float(rows[-1][0]), float(rows[-1][1]), float(rows[len(rows) // 2][1])


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
res = []
print(f"{'P (W)':>6} {'p (mTorr)':>9} {'Te (eV)':>22} {'n_e (1e16 m^-3)':>28} {'change, 2nd half':>22}")
for run in sorted(glob.glob(os.path.join(runs, "P*_p*"))):
    P, mT = (float(v) for v in re.fullmatch(r"P([0-9.]+)_p([0-9.]+)", os.path.basename(run)).groups())
    N = mT * 1e-3 * 133.3224 / (KB * Tg)
    t, Te, Te_mid = probe(run, "Te")
    _, ne, ne_mid = probe(run, "N_electron")
    TeL, neL = global_model(P, N)
    drift = max(abs(Te / Te_mid - 1), abs(ne / ne_mid - 1))
    res.append((P, mT, Te, ne, TeL, neL, drift))
    print(f"{P:6.0f} {mT:9.0f} {Te*KB/E:8.4f} ({TeL*KB/E:7.4f}, {Te/TeL-1:+.1e}) "
          f"{ne/1e16:10.5f} ({neL/1e16:9.5f}, {ne/neL-1:+.1e}) {drift:14.1e}")

ok = len(res) == 5 and all(abs(r[2] / r[4] - 1) < 1e-3 and abs(r[3] / r[5] - 1) < 1e-3 and r[6] < 1e-4 for r in res)

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    pw = sorted(r for r in res if r[1] == 10)
    PP = [10**(2 + i / 50) for i in range(51)]
    N10 = 10e-3 * 133.3224 / (KB * Tg)
    a1.loglog(PP, [global_model(p, N10)[1] for p in PP], "k-", lw=0.8, label="global model")
    a1.loglog([r[0] for r in pw], [r[3] for r in pw], "o", mfc="none", label="SOMAFOAM")
    a1.set_xlabel("absorbed power (W)")
    a1.set_ylabel("$n_\\mathrm{e}$ (m$^{-3}$)")
    a1.set_title("(a) 10 mTorr", fontsize=8)
    a1.legend(frameon=False, fontsize=7)
    pr = sorted((r for r in res if r[0] == 300), key=lambda r: r[1])
    mm = [10**(math.log10(2) + i * (math.log10(50) - math.log10(2)) / 50) for i in range(51)]
    a2.semilogx(mm, [global_model(300, m * 1e-3 * 133.3224 / (KB * Tg))[0] * KB / E for m in mm], "k-", lw=0.8)
    a2.semilogx([r[1] for r in pr], [r[2] * KB / E for r in pr], "o", mfc="none")
    a2.set_xlabel("pressure (mTorr)")
    a2.set_ylabel("$k_\\mathrm{B}T_\\mathrm{e}/e$ (V)")
    a2.set_title("(b) 300 W", fontsize=8)
    from matplotlib.ticker import NullFormatter, NullLocator
    for ax in (a1, a2):
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
    a1.yaxis.set_minor_formatter(NullFormatter())
    a1.set_xticks([100, 300, 1000])
    a1.set_xticklabels(["100", "300", "1000"])
    a2.set_xticks([3, 10, 30])
    a2.set_xticklabels(["3", "10", "30"])
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
