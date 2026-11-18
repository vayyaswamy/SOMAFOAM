"""Write the meshes and the boundary conditions of the coaxialDielectric2D case.

    python3 makeCase.py <case> <cells along the arc> <radial cells, dielectric> <radial cells, gas>

A quarter of a coaxial capacitor in the x-y plane, one cell thick: the
electrode at r = a = 1 mm, a dielectric (eps_r = 5, region dielectric_1) to
the curved interface at r = b = 2 mm, and gas to the grounded wall at
r = c = 5 mm. The cuts at x = 0 and y = 0 are symmetry planes. The meshes of
the two regions match on the interface.
"""
import math
import re
import sys

case, nt, nd, ng = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
a, b, c, h = 1e-3, 2e-3, 5e-3, 1e-4
V0 = 10.0
s45 = math.sqrt(0.5)

head = open(f"{case}/system/controlDict").read().split("FoamFile")[0]


def foamfile(cls, obj):
    return head + "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       %s;\n    object      %s;\n}\n\n" % (cls, obj)


def annulus(r1, r2, nr):
    v = [(r1, 0, 0), (r2, 0, 0), (0, r2, 0), (0, r1, 0), (r1, 0, h), (r2, 0, h), (0, r2, h), (0, r1, h)]
    s = "convertToMeters 1;\n\nvertices\n(\n" + "".join("    (%.12g %.12g %.12g)\n" % p for p in v) + ");\n\n"
    s += "blocks\n(\n    hex (0 1 2 3 4 5 6 7) (%d %d 1) simpleGrading (1 1 1)\n);\n\n" % (nr, nt)
    s += "edges\n(\n"
    for (i, j, r, z) in ((0, 3, r1, 0), (1, 2, r2, 0), (4, 7, r1, h), (5, 6, r2, h)):
        s += "    arc %d %d (%.12g %.12g %.12g)\n" % (i, j, r * s45, r * s45, z)
    s += ");\n\nboundary\n(\n"
    return s


def patch(name, body, faces):
    return "    %s\n    {\n%s        faces ( %s );\n    }\n" % (name, body, " ".join(faces))


common = [
    ("cutY", "        type symmetryPlane;\n", ["(0 1 5 4)"]),
    ("cutX", "        type symmetryPlane;\n", ["(3 7 6 2)"]),
    ("frontAndBack", "        type empty;\n", ["(0 3 2 1)", "(4 5 6 7)"]),
]
couple_master = """        type            regionCouple;
        shadowRegion    dielectric_1;
        shadowPatch     plasma_1;
        zone            plasma_1_Zone;
        attached        off;
        master          on;
        isWall          on;
        bridgeOverlap   off;
"""
couple_slave = """        type            regionCouple;
        shadowRegion    region0;
        shadowPatch     dielectric_1;
        zone            dielectric_1_Zone;
        attached        off;
        master          off;
        isWall          on;
        bridgeOverlap   off;
"""

# Gas (plasma region): b to c
s = foamfile("dictionary", "blockMeshDict") + annulus(b, c, ng)
s += patch("dielectric_1", couple_master, ["(0 4 7 3)"])
s += patch("ground", "        type patch;\n", ["(1 2 6 5)"])
for p in common:
    s += patch(*p)
s += ");\n\nmergePatchPairs ();\n"
open(f"{case}/constant/polyMesh/blockMeshDict", "w").write(s)

# Dielectric: a to b
s = foamfile("dictionary", "blockMeshDict") + annulus(a, b, nd)
s += patch("electrode", "        type patch;\n", ["(0 4 7 3)"])
s += patch("plasma_1", couple_slave, ["(1 2 6 5)"])
for p in common:
    s += patch(*p)
s += ");\n\nmergePatchPairs ();\n"
open(f"{case}/constant/dielectric_1/polyMesh/blockMeshDict", "w").write(s)


def entry(text, name):
    m = re.search(r"\n(\s*)%s\s*\n\s*\{" % re.escape(name), text)
    i = text.index("{", m.end() - 1)
    depth, j = 0, i
    while True:
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                break
        j += 1
    return text[i + 1:j]


def rewrite(path, entries, internal=None):
    t = open(path).read()
    if internal is not None:
        t = re.sub(r"internalField\s+[^;]+;", "internalField   %s;" % internal, t)
    k = t.index("boundaryField")
    body = "boundaryField\n{\n"
    for name, e in entries:
        body += "    %s\n    {%s    }\n" % (name, e)
    body += "}\n"
    open(path, "w").write(t[:k] + body)


sym = [("cutY", "\n        type symmetryPlane;\n"), ("cutX", "\n        type symmetryPlane;\n"),
       ("frontAndBack", "\n        type empty;\n")]

# Plasma region fields: the interface entries of the ArgonDBD example; the
# grounded wall takes the same wall conditions, the potential 0 and a
# zero-gradient permittivity
for f in ("electron", "Arp1", "Arm", "Ar", "p", "T", "Te", "Tion", "epsilon", "Phi"):
    path = f"{case}/0/{f}"
    t = open(path).read()
    e1 = entry(t, "dielectric_1")
    if f == "Phi":
        e1 = re.sub(r"surfaceCharge\s+\w+;", "surfaceCharge   off;", e1)
        eg = "\n        type            fixedValue;\n        value           uniform 0;\n"
    elif f == "epsilon":
        eg = "\n        type            zeroGradient;\n"
    else:
        eg = e1
    internal = "uniform 1e4" if f in ("electron", "Arp1", "Arm") else None
    rewrite(path, [("dielectric_1", e1), ("ground", eg)] + sym, internal)

# Dielectric fields
path = f"{case}/0/dielectric_1/Phi"
t = open(path).read()
rewrite(path, [("electrode", "\n        type            fixedValue;\n        value           uniform %g;\n" % V0),
               ("plasma_1", entry(t, "plasma_1"))] + sym, "uniform 0")
path = f"{case}/0/dielectric_1/epsilon"
t = open(path).read()
rewrite(path, [("electrode", "\n        type            zeroGradient;\n"),
               ("plasma_1", entry(t, "plasma_1"))] + sym)

# One dielectric region; a few steps of 0.1 ns
t = open(f"{case}/constant/regionProperties").read()
t = re.sub(r"dielectricRegionNames\s*\([^)]*\);", "dielectricRegionNames\n(\n    dielectric_1\n);", t)
open(f"{case}/constant/regionProperties", "w").write(t)
t = open(f"{case}/system/controlDict").read()
for k, v in (("endTime", "5e-10"), ("writeInterval", "5e-10")):
    t = re.sub(r"\n%s\s+[^;]+;" % k, "\n%-15s %s;" % (k, v), t)
open(f"{case}/system/controlDict", "w").write(t)
