"""Build an axisymmetric argon ICP case for somaFoam (electroMagnetics dischargeType ICP).

Usage:
  python3 make_case.py TEMPLATE DBD_TEMPLATE OUTDIR [--window wall|dielectric] [--refine K]
                       [--nproc N] [--power W] [--pressure mTorr] [--vcoil V] [--verify]

TEMPLATE is the 1 Torr Scharfetter-Gummel argon example (plasma-only), DBD_TEMPLATE the
ArgonDBD example (for the dielectric-region settings). The planar-coil geometry matches the
preliminary OpenFOAM-10 and finite-difference solutions: plasma 0 < r < R = 10 cm,
0 < z < L = 10 cm, quartz window L < z < 11.2 cm, five-turn coil of square 4 mm conductors
at z = 11.85 cm, metal top at 16 cm. With --window dielectric the window is a dielectric
region (epsilon_r = 3.8) coupled to the plasma, and the part of its top surface under the
coil carries the RF coil voltage, which supplies the capacitive (E-mode) coupling.
With --verify the ICP module imposes the analytic density n0 J0(2.405 r/R) sin(pi z/L) of
the preliminary solutions, solves the field once, prints the power per A^2 and exits.
"""
import sys, os, re, shutil, math, argparse

ap = argparse.ArgumentParser()
ap.add_argument("tpl"); ap.add_argument("dbd"); ap.add_argument("out")
ap.add_argument("--window", default="wall", choices=["wall", "dielectric"])
ap.add_argument("--refine", type=int, default=1)
ap.add_argument("--nproc", type=int, default=1)
ap.add_argument("--power", type=float, default=30.0)
ap.add_argument("--pressure", type=float, default=50.0)
ap.add_argument("--vcoil", type=float, default=300.0)
ap.add_argument("--freq", type=float, default=13.56e6)
ap.add_argument("--verify", action="store_true")
ap.add_argument("--rates", default="maxwellian", choices=["maxwellian", "template"],
                help="electron-impact rate tables: Maxwellian fits (ICP bulk) or the template's BOLSIG+ tables")
ap.add_argument("--meta-wall", default="quench", choices=["quench", "zeroGradient"],
                help="metastable wall condition: quench (n = 0, de-excitation probability near one) or the template's zeroGradient")
ap.add_argument("--meta-accel", type=float, default=1.0, help="time acceleration of the metastable equation")
ap.add_argument("--end-time", type=float, default=2e-4)
ap.add_argument("--dt-max", default="1e-10", help="largest time step [s]")
ap.add_argument("--write-interval", type=float, default=5e-6)
A = ap.parse_args()
k = A.refine
R, L, ZWIN, ZTOP = 0.10, 0.10, 0.112, 0.16
diel = A.window == "dielectric"

hdr = lambda cls, obj, loc: ("FoamFile\n{\n    version 2.0;\n    format ascii;\n    class %s;\n"
                             "    location \"%s\";\n    object %s;\n}\n" % (cls, loc, obj))
# radial blocks: axis, under the coil, wall layer; axial blocks: wall layer, bulk, window layer
rb, rn, rg = [0, 0.015, 0.085, R], [5*k, 21*k, 14*k], [1, 1, 1/50.]
zb_p, zn_p, zg_p = [0, 0.015, 0.085, L], [12*k, 20*k, 12*k], [50, 1, 1/50.]
zb_w, zn_w, zg_w = [L, ZWIN], [4*k], [1]
zb_a, zn_a, zg_a = [ZWIN, 0.124, ZTOP], [6*k, 8*k], [1, 1]
a = math.radians(2.5)

def bmd(zb, zn, zg, patchfn, loc):
    """Wedge blockMeshDict; patchfn(i, j, side) returns the patch of a block face."""
    V, idx = [], {}
    def v(i, j, s):
        key = (i, j, s)
        if key not in idx:
            r, z = rb[i], zb[j]
            idx[key] = len(V); V.append((r*math.cos(a), s*r*math.sin(a), z))
        return idx[key]
    blocks, P = [], {}
    for j in range(len(zn)):
        for i in range(len(rn)):
            h = [v(i,j,-1), v(i+1,j,-1), v(i+1,j,1), v(i,j,1), v(i,j+1,-1), v(i+1,j+1,-1), v(i+1,j+1,1), v(i,j+1,1)]
            blocks.append("hex (%s) (%d 1 %d) simpleGrading (%g 1 %g)" % (" ".join(map(str, h)), rn[i], zn[j], rg[i], zg[j]))
            P.setdefault("front", []).append((h[0],h[1],h[5],h[4])); P.setdefault("back", []).append((h[3],h[7],h[6],h[2]))
            if i == 0: P.setdefault("axis", []).append((h[0],h[4],h[7],h[3]))
            if i == len(rn)-1: P.setdefault(patchfn(i, j, "outer"), []).append((h[1],h[2],h[6],h[5]))
            if j == 0: P.setdefault(patchfn(i, j, "bottom"), []).append((h[0],h[3],h[2],h[1]))
            if j == len(zn)-1: P.setdefault(patchfn(i, j, "top"), []).append((h[4],h[5],h[6],h[7]))
    f = lambda fl: " ".join("(%d %d %d %d)" % t for t in fl)
    def ptype(n):
        if n in ("front", "back"): return "wedge"
        if n == "axis": return "empty"
        if n == "dielectric_1":
            return ("regionCouple; shadowRegion dielectric_1; shadowPatch plasma_1; zone plasma_1_Zone; "
                    "attached off; master on; isWall on; bridgeOverlap off")
        if n == "plasma_1":
            return ("regionCouple; shadowRegion region0; shadowPatch dielectric_1; zone dielectric_1_Zone; "
                    "attached off; master off; isWall on; bridgeOverlap off")
        return "patch"
    body = "\n".join("    %s { type %s; faces (%s); }" % (n, ptype(n), f(fl)) for n, fl in P.items())
    return (hdr("dictionary", "blockMeshDict", loc) + "convertToMeters 1;\nvertices\n(\n" +
            "\n".join("    (%.12g %.12g %.12g)" % p for p in V) + "\n);\nblocks\n(\n" +
            "\n".join("    " + b for b in blocks) + "\n);\nedges ();\nboundary\n(\n" + body + "\n);\nmergePatchPairs ();\n")

out = A.out
if os.path.exists(out): shutil.rmtree(out)
shutil.copytree(A.tpl, out)
shutil.rmtree(os.path.join(out, "constant/polyMesh"))
def wfile(rel, txt):
    p = os.path.join(out, rel); os.makedirs(os.path.dirname(p), exist_ok=True); open(p, "w").write(txt)

# plasma mesh: window face is a wall (grounded) or the coupled dielectric interface
wfile("constant/polyMesh/blockMeshDict", bmd(zb_p, zn_p, zg_p,
      lambda i, j, s: ("dielectric_1" if (diel and s == "top") else "wall"), "constant/polyMesh"))
# em mesh: plasma blocks first (identical cell order), then window and air
wfile("constant/em/polyMesh/blockMeshDict", bmd(zb_p + zb_w[1:] + zb_a[1:], zn_p + zn_w + zn_a, zg_p + zg_w + zg_a,
      lambda i, j, s: "wall", "constant/em/polyMesh"))
if diel:
    wfile("constant/dielectric_1/polyMesh/blockMeshDict", bmd(zb_w, zn_w, zg_w,
          lambda i, j, s: {"bottom": "plasma_1", "outer": "wall"}.get(s, "coil" if i == 1 else "top"),
          "constant/dielectric_1/polyMesh"))

def entry(txt, name):
    m = re.search(r"\n\s*%s\s*\n\s*\{" % name, txt); i = txt.index("{", m.start()); d = 0
    for q in range(i, len(txt)):
        d += txt[q] == "{"; d -= txt[q] == "}"
        if d == 0: return txt[i:q+1]
set_ic = {"p": A.pressure*0.133322, "electron": 1e16, "Arp1": 1e16, "Arm": 1e16, "Te": 23000}
for fn in os.listdir(os.path.join(out, "0")):
    p = os.path.join(out, "0", fn)
    if os.path.isdir(p): continue
    t = open(p).read()
    if "boundaryField" not in t: continue
    g = entry(t, "ground").replace("\n", "\n    ")
    win = g
    if diel and fn == "Phi":
        win = "{ type coupledPotential; surfaceCharge on; remoteField Phi; value uniform 0; }"
    t = t[:t.index("boundaryField")] + ("boundaryField\n{\n    wall\n    %s\n%s    front { type wedge; }\n    back { type wedge; }\n"
        "    axis { type empty; }\n}\n" % (g, ("    dielectric_1\n    %s\n" % win) if diel else ""))
    if fn == "Arm" and A.meta_wall == "quench":
        # at 50 mTorr the extrapolation length D/(v_th/4) of a wall with de-excitation probability
        # near one is about 1 mm, so a Dirichlet condition represents wall quenching
        t = re.sub(r"(\n    wall\n    )\{[^}]*\}", r"\1{ type fixedValue; value uniform 1e6; }", t)
    if fn in set_ic: t = re.sub(r"internalField\s+uniform\s+\S+;", "internalField   uniform %g;" % set_ic[fn], t)
    open(p, "w").write(t)

if diel:
    eps0 = 8.8541878128e-12
    wfile("0/epsilon", hdr("volScalarField", "epsilon", "0") + "dimensions [0 0 0 1 0 0 0];\ninternalField uniform %g;\n"
          "boundaryField\n{\n    wall { type zeroGradient; }\n    dielectric_1 { type plasmaDielectricEpsilon; remoteField epsilon; value uniform %g; }\n"
          "    front { type wedge; }\n    back { type wedge; }\n    axis { type empty; }\n}\n" % (eps0, eps0))
    wfile("0/dielectric_1/epsilon", hdr("volScalarField", "epsilon", "0/dielectric_1") + "dimensions [0 0 0 1 0 0 0];\n"
          "internalField uniform %g;\nboundaryField\n{\n    plasma_1 { type plasmaDielectricEpsilonSlave; remoteField epsilon; value uniform %g; }\n"
          "    coil { type zeroGradient; }\n    top { type zeroGradient; }\n    wall { type zeroGradient; }\n"
          "    front { type wedge; }\n    back { type wedge; }\n    axis { type empty; }\n}\n" % (3.8*eps0, 3.8*eps0))
    wfile("0/dielectric_1/Phi", hdr("volScalarField", "Phi", "0/dielectric_1") + "dimensions [0 0 0 1 0 0 0];\n"
          "internalField uniform 0;\nboundaryField\n{\n    plasma_1 { type coupledPotential; surfaceCharge off; remoteField Phi; value uniform 0; }\n"
          "    coil\n    {\n        type plasmaPotential;\n        value uniform 0;\n        model sinFrequencyModulated;\n"
          "        amplitude uniform %g;\n        frequency %g;\n        bias 0;\n        R 0;\n    }\n"
          "    top { type zeroGradient; }\n    wall { type fixedValue; value uniform 0; }\n"
          "    front { type wedge; }\n    back { type wedge; }\n    axis { type empty; }\n}\n" % (A.vcoil, A.freq))
    os.makedirs(os.path.join(out, "system/dielectric_1"), exist_ok=True)
    for f in ("fvSchemes", "fvSolution"):
        shutil.copy(os.path.join(A.dbd, "system/dielectric_1", f), os.path.join(out, "system/dielectric_1", f))
    wfile("constant/regionProperties", hdr("dictionary", "regionProperties", "constant") + "dielectricRegionNames\n(\n    dielectric_1\n);\n")
    p = os.path.join(out, "system/fvSolution"); t = open(p).read()
    for key in ("Phi_Phi_Phi", "Phi_Phi"):
        if not re.search(r"\n\s*%s\s*\n" % key, t):
            i = t.index("solvers"); i = t.index("{", i) + 1
            t = t[:i] + "\n    %s\n    {\n        solver BiCGStab;\n        preconditioner Cholesky;\n        relTol 0;\n        minIter 1;\n    }\n" % key + t[i:]
    open(p, "w").write(t)
    wfile("dielectric_1.setSet", "faceSet dielectric_1_Zone new patchToFace plasma_1;\nquit;\n")
    wfile("plasma_1.setSet", "faceSet plasma_1_Zone new patchToFace dielectric_1;\nquit;\n")

# em-region fields and numerics
bf = "boundaryField\n{\n    wall { type fixedValue; value uniform %s; }\n    front { type wedge; }\n    back { type wedge; }\n    axis { type empty; }\n}\n"
for n in ("AR", "AI"):
    wfile("0/em/" + n, hdr("volScalarField", n, "0/em") + "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\n" + bf % "0")
wfile("system/em/fvSchemes", hdr("dictionary", "fvSchemes", "system/em") +
      "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\ndivSchemes { default none; }\n"
      "laplacianSchemes { default Gauss linear corrected; }\ninterpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
wfile("system/em/fvSolution", hdr("dictionary", "fvSolution", "system/em") +
      "solvers\n{\n    A { solver BiCGStab; preconditioner Cholesky; tolerance 1e-10; relTol 0; minIter 1; maxIter 20000; }\n}\n")

# decomposition: radial strips with equal cell counts give the same radial cuts in every
# region, so each processor holds the em and dielectric cells above its plasma cells
dec = lambda zones: (hdr("dictionary", "decomposeParDict", "system") + "numberOfSubdomains %d;\nmethod simple;\n"
                     "simpleCoeffs { n (%d 1 1); delta 0.001; }\n%sdistributed no;\nroots ();\n"
                     % (A.nproc, A.nproc, ("globalFaceZones (%s);\n" % zones) if zones else ""))
wfile("system/decomposeParDict", dec("plasma_1_Zone" if diel else ""))
wfile("system/em/decomposeParDict", dec(""))
if diel: wfile("system/dielectric_1/decomposeParDict", dec("dielectric_1_Zone"))

# electroMagnetics: ICP option
p = os.path.join(out, "constant/electroMagnetics"); t = open(p).read()
if diel: t = re.sub(r"solutionDomain\s+\w+;", "solutionDomain\t\tplasmaDielectric;", t)
t += ("\ndischargeType   ICP;\nICPCoeffs\n{\n    frequency %g;\n    power %g;\n    wedgeAngle 5;\n    updateInterval 20;\n"
      "    updateTolerance 0.01;\n    maxUpdateInterval 10000;\n    powerRampTime 5e-6;\n    powerRelaxation 0.5;\n    coilR (0.02 0.035 0.05 0.065 0.08);\n    coilZ 0.1185;\n    coilW 0.004;\n%s}\n"
      % (A.freq, A.power, "    verification { n0 5.0766242253677146e+17; nuM 33573101.15215304; R 0.10; L 0.10; }\n" if A.verify else ""))
open(p, "w").write(t)

# plasmaProperties: a metastable given with diffusionModel constant (a D at the template pressure)
# is converted to diffusionModel pressureScaled with pD = D p_template, so that D = pD/p follows
# the gas pressure as in the CCP cases; the charged-species tables are density-reduced already.
# Optional time acceleration of the metastable equation.
p_tpl = float(re.search(r"internalField\s+uniform\s+(\S+);", open(os.path.join(A.tpl, "0/p")).read()).group(1))
p = os.path.join(out, "constant/plasmaProperties"); t = open(p).read()
i = t.index("\nArm"); j = t.index("}", i)
blk = t[i:j]
if re.search(r"diffusionModel\s+constant;", blk):
    D0 = float(re.search(r"\n\s*D\s+(\S+);", blk).group(1)); pD = D0*p_tpl
    blk = re.sub(r"diffusionModel\s+constant;", "diffusionModel pressureScaled;", blk)
    blk = re.sub(r"(\n\s*)D\s+\S+;", lambda m: m.group(1) + "pD %.6g;" % pD, blk)
else:
    pD = float(re.search(r"\n\s*pD\s+(\S+);", blk).group(1)); D0 = pD/p_tpl
Dn = pD/(A.pressure*0.133322)
blk = re.sub(r"\n\s*timeAcceleration\s+\S+;", "", blk) + "\ttimeAcceleration %g;\n" % A.meta_accel
t = t[:i] + blk + t[j:]; open(p, "w").write(t)
print("Arm D %.4g -> %.4g m2/s (template p %.4g Pa), timeAcceleration %g" % (D0, Dn, p_tpl, A.meta_accel))

# Rate tables. At ICP densities electron-electron collisions bring the distribution close to a
# Maxwellian, so by default the ground-state ionization and excitation use the Maxwellian fits of
# Lieberman and Lichtenberg, and stepwise ionization of Ar* (threshold 4.21 eV) uses the Maxwellian
# average of the Lotz cross section 4.5e-18 m^2 eV^2 ln(E/I)/(E I). The total excitation fit is used
# as the metastable source, which is an upper bound. The template's BOLSIG+ tables remain available
# with --rates template as a sensitivity case. Tables hold ln(k) with k in m^3/(kmol s) against Te in K.
if A.rates == "maxwellian":
    import numpy as _np
    NAk = 6.02214076e26; eVK = 11604.518
    tpl_tab = open(os.path.join(A.tpl, "constant/reaction_2")).read()
    TK = _np.array([float(a) for a, b in re.findall(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", tpl_tab)])
    Te = TK/eVK
    def lotz(T, I=4.21, a=4.5e-18):
        E = _np.linspace(I, I + 60*max(T, 0.05), 4000)
        sig = a*_np.log(E/I)/(E*I)
        v = _np.sqrt(2*E*1.602176634e-19/9.1093837e-31)
        f = 2*_np.sqrt(E/_np.pi)*T**-1.5*_np.exp(-E/T)  # Maxwellian energy distribution, 1/eV
        return (_np.trapezoid if hasattr(_np, 'trapezoid') else _np.trapz)(sig*v*f, E)
    K = {"reaction_2": 2.34e-14*Te**0.59*_np.exp(-17.44/Te),
         "reaction_3": 2.48e-14*Te**0.33*_np.exp(-12.78/Te),
         "reaction_met": _np.array([lotz(t) for t in Te])}
    for fn, kt in K.items():
        rows = "\n".join("(%.10g %.10g)" % (x, _np.log(max(kk, 1e-40)*NAk)) for x, kk in zip(TK, kt))
        open(os.path.join(out, "constant", fn), "w").write("(\n" + rows + "\n)\n")
    for t in (2.0, 3.0, 3.6):
        i = _np.argmin(abs(Te - t))
        print("Maxwellian rates at %.2f eV: iz %.2e, ex %.2e, stepwise iz %.2e m3/s"
              % (Te[i], K["reaction_2"][i], K["reaction_3"][i], K["reaction_met"][i]))

# controlDict: adaptive step, wall-current output, long run (RF-resolving step with the coil voltage)
p = os.path.join(out, "system/controlDict"); t = open(p).read()
dtmax = A.dt_max  # 1e-10 s, the step of the 1 Torr template, by default
for key, val in [("endTime", "%g" % A.end_time), ("deltaT", "1e-10"), ("deltaTMax", dtmax), ("deltaTMin", "1e-12"),
                 ("writeInterval", "%g" % A.write_interval), ("MaxCo", "10.0")]:
    t = re.sub(r"\n%s\s+[^;]*;" % key, "\n%s    %s;" % (key, val), t)
t = re.sub(r"patches\s+\(electrode ground\);", "patches             (wall);", t)
open(p, "w").write(t)

# mesh script
steps = []
if diel:
    steps += ["blockMesh -region dielectric_1 > log.blockMesh.dielectric_1 2>&1",
              "setSet -region dielectric_1 -batch dielectric_1.setSet -noVTK > log.setSet.dielectric_1 2>&1",
              "setsToZones -region dielectric_1 -noFlipMap > log.setsToZones.dielectric_1 2>&1"]
steps += ["blockMesh > log.blockMesh 2>&1"]
if diel:
    steps += ["setSet -batch plasma_1.setSet -noVTK > log.setSet 2>&1", "setsToZones -noFlipMap > log.setsToZones 2>&1"]
steps += ["blockMesh -region em > log.blockMesh.em 2>&1"]
if A.nproc > 1:
    steps += ["decomposePar > log.decomposePar 2>&1", "decomposePar -region em > log.decomposePar.em 2>&1"]
    if diel: steps += ["decomposePar -region dielectric_1 > log.decomposePar.dielectric_1 2>&1"]
wfile("Allmesh", "#!/bin/bash\nset -e\ncd ${0%/*}\n" + "\n".join(s + ' || { echo "failed: ' + s.split()[0] + '"; exit 1; }' for s in steps) + "\necho mesh-ok\n")
os.chmod(os.path.join(out, "Allmesh"), 0o755)
print("case written to", out, "| window", A.window, "| refine", k, "| nproc", A.nproc, "| verify", A.verify)
