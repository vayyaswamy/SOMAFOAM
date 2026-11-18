"""Write the wedge mesh and the solenoid of the icpSkinEffect case.

    python3 makeCase.py <case> <radial cells> [wedge angle, degrees, default 1]

Axisymmetric wedge (1 degree by default) about the z axis, from r = 1e-5 m (axis
patch, A = 0) to the conducting wall at r = 0.1 m in two blocks that meet at
the edge of the plasma, r = 0.05 m, and 50 mm long with
zero-gradient ends, so that the solenoid of 20 contiguous square turns of
2.5 mm at r = 0.07 m is infinitely long.
"""
import math
import os
import sys

case, nr = sys.argv[1], int(sys.argv[2])
ang = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
r0, rw, H, nz = 1e-5, 0.1, 0.05, 20
rc, w, nturn = 0.07, 2.5e-3, 20
c, s = math.cos(math.radians(ang / 2)), math.sin(math.radians(ang / 2))

rp = 0.05
# two radial blocks, so that the edge of the plasma at r = 0.05 m is a face
v = []
for z in (0, H):
    for r in (r0, rp, rw):
        v.append((r * c, -r * s, z))
    for r in (r0, rp, rw):
        v.append((r * c, r * s, z))
# z = 0: 0 1 2 (-y), 3 4 5 (+y); z = H: 6 7 8 (-y), 9 10 11 (+y)

hdr = open(f"{case}/system/controlDict").read().split("FoamFile")[0]
os.makedirs(f"{case}/constant/polyMesh", exist_ok=True)
with open(f"{case}/constant/polyMesh/blockMeshDict", "w") as f:
    f.write(hdr + "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      blockMeshDict;\n}\n\n")
    f.write("convertToMeters 1;\n\nvertices\n(\n")
    for p in v:
        f.write("    (%.12g %.12g %.12g)\n" % p)
    f.write(");\n\nblocks\n(\n")
    f.write("    hex (0 1 4 3 6 7 10 9) (%d 1 %d) simpleGrading (1 1 1)\n" % (nr // 2, nz))
    f.write("    hex (1 2 5 4 7 8 11 10) (%d 1 %d) simpleGrading (1 1 1)\n" % (nr // 2, nz))
    f.write(");\n\nedges ();\n\n")
    f.write("""patches
(
    patch axis   ( (0 6 9 3) )
    patch wall   ( (2 5 11 8) )
    patch ends   ( (0 3 4 1) (1 4 5 2) (6 7 10 9) (7 8 11 10) )
    wedge front  ( (0 1 7 6) (1 2 8 7) )
    wedge back   ( (3 9 10 4) (4 10 11 5) )
);

mergePatchPairs ();
""")

turns = " ".join("(%.6g %.6g)" % (rc, w / 2 + k * w) for k in range(nturn))
with open(f"{case}/system/icpFieldDict", "w") as f:
    f.write(hdr + "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      icpFieldDict;\n}\n\n")
    f.write(f"""frequency   13.56e6;
n0          1e17;
nuM         1e7;
profile     uniform;
R           0.05;
coilW       {w};
coilTurns   ({turns});
""")
