"""Vector potential of the icpSkinEffect runs against the exact solution.

    python3 verify.py runs [--plot figure.pdf]

A uniform plasma of radius Rp (n0 = 1e17 m^-3, nu_m = 1e7 1/s, 13.56 MHz)
inside an infinitely long solenoid of contiguous square turns of side w at
radius Rc, 1 A per turn, and a conducting wall at Rw where A = 0. With
kappa^2 = i w mu0 sigma and sigma = n0 e^2/(m (nu_m + i w)), the axial field
is B1 I0(kappa r)/I0(kappa Rp) in the plasma, B1 between the plasma and the
winding, falls linearly by mu0/w across the winding, and is B1 - mu0/w
outside it; B1 follows from the zero total flux inside the wall, and
A = flux(r)/(2 pi r). The error is the largest |A - A_exact| over the cells
relative to the largest |A_exact|.
"""
import cmath
import glob
import math
import os
import re
import sys

f, n0, nu = 13.56e6, 1e17, 1e7
Rp, Rc, w, Rw = 0.05, 0.07, 2.5e-3, 0.1
qe, me, mu0 = 1.602176634e-19, 9.1093837e-31, 4e-7 * math.pi
om = 2 * math.pi * f
sigma = n0 * qe * qe / (me * (nu + 1j * om))
kappa = cmath.sqrt(1j * om * mu0 * sigma)
a, b = Rc - w / 2, Rc + w / 2
j = 1.0 / w**2                      # current density of a turn [A/m2]


def besselI(nu_, x, terms=80):
    s, t = 0, (x / 2)**nu_ / math.factorial(nu_)
    for m in range(terms):
        s += t
        t *= (x / 2)**2 / ((m + 1) * (m + 1 + nu_))
    return s


I0p = besselI(0, kappa * Rp)


def flux(r, B1):
    """Axial flux inside radius r for the field B1 between plasma and winding."""
    Fp = 2 * math.pi * B1 * Rp * besselI(1, kappa * Rp) / (kappa * I0p)
    if r <= Rp:
        return 2 * math.pi * B1 * r * besselI(1, kappa * r) / (kappa * I0p)
    if r <= a:
        return Fp + B1 * math.pi * (r * r - Rp * Rp)
    Fa = Fp + B1 * math.pi * (a * a - Rp * Rp)
    if r <= b:
        return Fa + 2 * math.pi * (B1 * (r * r - a * a) / 2 - mu0 * j * ((r**3 - a**3) / 3 - a * (r * r - a * a) / 2))
    Fb = Fa + 2 * math.pi * (B1 * (b * b - a * a) / 2 - mu0 * j * ((b**3 - a**3) / 3 - a * (b * b - a * a) / 2))
    B2 = B1 - mu0 * j * w
    return Fb + B2 * math.pi * (r * r - b * b)


# flux(Rw) is linear in B1: flux = B1 f1 + f0 = 0
f0 = flux(Rw, 0)
f1 = flux(Rw, 1) - f0
B1 = -f0 / f1


def Aexact(r):
    return flux(r, B1) / (2 * math.pi * r)


def cells(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    return [float(v) for v in s[m.end():].split(")")[0].split()][:int(m.group(1))]


runs = sys.argv[1]
plot = sys.argv[sys.argv.index("--plot") + 1] if "--plot" in sys.argv else None
results = []
for d in sorted(glob.glob(os.path.join(runs, "nr_*")), key=lambda p: int(p.split("_")[-1])):
    nr = int(d.split("_")[-1])
    AR, AI = cells(os.path.join(d, "1", "AR")), cells(os.path.join(d, "1", "AI"))
    nz = len(AR) // nr
    r0 = 1e-5
    # radial faces of the two blocks, and the centroid radius of each cell
    faces = [r0 + i * (Rp - r0) / (nr // 2) for i in range(nr // 2)] + [Rp + i * (Rw - Rp) / (nr // 2) for i in range(nr // 2 + 1)]
    rcell = [2 / 3 * (b_**3 - a_**3) / (b_**2 - a_**2) for a_, b_ in zip(faces, faces[1:])]
    num = []
    zs = []
    for i in range(nr):
        # blockMesh numbers the cells block by block, radius fastest
        h = nr // 2
        idx = [(k * h + i) if i < h else (h * nz + k * h + i - h) for k in range(nz)]
        col = [complex(AR[m], AI[m]) for m in idx]
        num.append(sum(col) / nz)
        zs.append(max(abs(c - num[-1]) for c in col))
    ex = [Aexact(r) for r in rcell]
    scale = max(abs(e) for e in ex)
    err = max(abs(n_ - e) for n_, e in zip(num, ex)) / scale
    zdev = max(zs) / scale
    log = open(os.path.join(d, "log.icpFieldTest")).read()
    P = float(re.search(r"power per \(1 A\)\^2 in the wedge (\S+) W", log).group(1))
    results.append((nr, err, zdev, P, rcell, num, ex))
    print(f"{nr:4d} radial cells: max |A - A_exact|/max|A_exact| = {err:.3e}, "
          f"variation along z {zdev:.1e}, power per A^2 (full circle) {P*360/5:.6e} W")

# Power per square ampere of the exact solution, 1/2 Re(sigma) w^2 int |A|^2 2 pi r dr
nint = 20000
Pex = 0
for i in range(nint):
    r = (i + 0.5) * Rp / nint
    Pex += 0.5 * sigma.real * om**2 * abs(Aexact(r))**2 * 2 * math.pi * r * Rp / nint
Pex *= 0.05                          # length of the solenoid section [m]
print(f"exact power per A^2 in 50 mm of length {Pex:.6e} W; skin depth "
      f"{1/abs(kappa.real):.4e} m, |kappa Rp| = {abs(kappa*Rp):.3f}")

ok = len(results) >= 2
for (n1, e1, *_), (n2, e2, *_) in zip(results, results[1:]):
    p = math.log(e1 / e2) / math.log(n2 / n1)
    print(f"order between {n1} and {n2} cells: {p:.2f}")
    ok &= p > 1.8
for r in results:
    print(f"{r[0]:4d} cells: power error {abs(r[3]*72/Pex - 1):.2e}")
ok &= results[-1][1] < 1e-3 and abs(results[-1][3] * 72 / Pex - 1) < 1e-3

if plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 2.8))
    nr, err, zdev, P, rc, num, ex = results[-1]
    rr = [i * Rw / 400 + 1e-6 for i in range(401)]
    exf = [Aexact(r) for r in rr]
    a1.plot([r * 100 for r in rr], [e.real * 1e6 for e in exf], "k-", lw=0.8, label="exact")
    a1.plot([r * 100 for r in rr], [e.imag * 1e6 for e in exf], "k--", lw=0.8)
    a1.plot([r * 100 for r in rc[::25]], [n_.real * 1e6 for n_ in num[::25]], "o", ms=3, mfc="none", label="Re A")
    a1.plot([r * 100 for r in rc[::25]], [n_.imag * 1e6 for n_ in num[::25]], "s", ms=3, mfc="none", label="Im A")
    a1.axvspan(0, Rp * 100, color="0.92")
    a1.axvspan(a * 100, b * 100, color="0.8")
    a1.set_xlabel("r (cm)")
    a1.set_ylabel("A per ampere-turn ($\\mu$Wb/m)")
    a1.legend(frameon=False, fontsize=7)
    a1.set_title("(a) vector potential, plasma and coil shaded", fontsize=8)
    a2.loglog([r[0] for r in results], [r[1] for r in results], "o-", ms=4, label="max error of A")
    a2.loglog([r[0] for r in results], [abs(r[3] * 72 / Pex - 1) for r in results], "s--", ms=4, label="error of the power")
    x0 = results[0][0]
    a2.loglog([r[0] for r in results], [results[0][1] * (x0 / r[0])**2 for r in results], "k:", label="$N^{-2}$")
    a2.set_xlabel("radial cells")
    a2.set_ylabel("relative error")
    a2.legend(frameon=False, fontsize=7)
    a2.set_title("(b) convergence", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot)
    print("wrote", plot)

print("PASS" if ok else "FAIL")
