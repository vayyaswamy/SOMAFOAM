"""Move the interior points of the unit-square mesh of a case.

    python3 distortMesh.py <case> <amplitude>

x and y are both displaced by a sin(2 pi x) sin(2 pi y). The boundary does
not move, and the displacement is smooth, so the cells become skewed and
non-orthogonal while the mesh stays a smooth map of the uniform one.
"""
import math
import re
import sys

case, a = sys.argv[1], float(sys.argv[2])
path = f"{case}/constant/polyMesh/points"
text = open(path).read()

m = re.search(r"\n(\d+)\s*\n\(", text)
n = int(m.group(1))
head, body = text[:m.end()], text[m.end():]
pts = re.findall(r"\(([^()]*)\)", body)[:n]
tail = body[body.index(pts[-1]) + len(pts[-1]) + 1:]

out = []
for p in pts:
    x, y, z = (float(v) for v in p.split())
    d = a * math.sin(2 * math.pi * x) * math.sin(2 * math.pi * y)
    out.append(f"({x + d:.15g} {y + d:.15g} {z:.15g})")

open(path, "w").write(head + "\n" + "\n".join(out) + tail)
