"""Compare a gas-flow verification case (examples/flow) with its analytic solution.

    python3 verify.py <channelStartup|advectionDiffusion> <case directory>
"""
import math, os, re, sys

L = 0.02
NX = 200


def latest(case):
    times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d) and float(d) > 0]
    return os.path.join(case, max(times, key=float))


def internal(path, vector=False):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", s)
    if not m:
        u = re.search(r"internalField\s+uniform\s+(\([^)]*\)|[0-9.eE+-]+)", s).group(1)
        v = [float(x) for x in u.strip("()").split()] if vector else float(u)
        return [v] * NX
    body = s[m.end():]
    if vector:
        vals = re.findall(r"\(([^()]*)\)", body)[:int(m.group(1))]
        return [[float(x) for x in v.split()] for v in vals]
    return [float(x) for x in body.split(")")[0].split()][:int(m.group(1))]


def report(name, x, num, exact, scale, tol):
    err = max(abs(a - b) for a, b in zip(num, exact)) / scale
    print(f"{name}: max |error| = {err:.3e} of the scale {scale:.4g} (tolerance {tol:g})")
    for i in (0, len(x) // 4, len(x) // 2):
        print(f"    x = {x[i]*1e3:6.3f} mm: computed {num[i]:.6g}, exact {exact[i]:.6g}")
    print("PASS" if err < tol else "FAIL")
    return err < tol


case, d = sys.argv[1], sys.argv[2]
t = latest(d)
x = [(i + 0.5) * L / NX for i in range(NX)]

p0, T0, W, R = 133.3224, 300.0, 40.0, 8314.47
rho = p0 * W / (R * T0)
core = [i for i, xi in enumerate(x) if 0.25 * L < xi < 0.75 * L]

if case == "channelStartup":
    f = 0.418
    U = internal(f"{t}/U", vector=True)[:NX]
    exact = [f * float(os.path.basename(t)) / rho] * NX
    ok = report("U_y, uniformly accelerated gas (core)", [x[i] for i in core], [U[i][1] for i in core],
                [exact[i] for i in core], exact[0], 1e-3)
elif case == "heatingStartup":
    Q, cp = 1e4, 520.3
    cv = cp - R / W
    T = internal(f"{t}/T")
    dT = Q * float(os.path.basename(t)) / (rho * cv)
    ok = report("T - 300 K, heating at constant volume (core)", [x[i] for i in core], [T[i] - T0 for i in core],
                [dT] * len(core), dT, 1e-2)
elif case == "advectionDiffusion":
    U0, D, N0, NL = 1.0, 4.7082e-3, 1e16, 1e10
    Pe = U0 * L / D
    N = internal(f"{t}/Arm")
    exact = [N0 + (NL - N0) * (math.exp(Pe * xi / L) - 1) / (math.exp(Pe) - 1) for xi in x]
    ok = report(f"N_m, advection-diffusion (Pe = {Pe:.3f})", x, N, exact, N0, 5e-3)
else:
    sys.exit(f"unknown case {case}")

sys.exit(0 if ok else 1)
