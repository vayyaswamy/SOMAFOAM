"""Mean energy and mobility of the boltzmannModelGas runs against exact.py.

    python3 verify.py runs

reads constant/boltzmann/transport.dat of each run and prints the relative
errors of the mean energy and of mu N at each reduced field, and the order
of the error in the number of energy cells.
"""
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact import exact

runs = sys.argv[1]
ref = {}
res = {}
for d in sorted(glob.glob(os.path.join(runs, "*_*"))):
    gas, n = re.fullmatch(r"(\w+)_(\d+)", os.path.basename(d)).groups()
    lines = open(os.path.join(d, "constant", "boltzmann", "transport.dat")).read().splitlines()
    head = [l for l in lines if l.startswith("#")]
    rows = [[float(x) for x in l.split()] for l in lines if l.strip() and not l.startswith("#")]
    for r in rows:
        en = r[0]
        if (gas, en) not in ref:
            ref[(gas, en)] = exact(gas, en)
        m, mu = ref[(gas, en)]
        res.setdefault((gas, en), []).append((int(n), r[1] / m - 1, r[3] / mu - 1))

ok = True
print(head[-1] if head else "")
for (gas, en), v in sorted(res.items()):
    v.sort()
    line = f"{gas:13s} {en:4g} Td: exact mean energy {ref[(gas, en)][0]:.5f} eV, mu N {ref[(gas, en)][1]:.5e};"
    for n, em, eu in v:
        line += f"  {n}: {em:+.1e} {eu:+.1e}"
    if len(v) > 2:
        p = math.log(abs(v[-2][1]) / abs(v[-1][1])) / math.log(v[-1][0] / v[-2][0]) if v[-1][1] and v[-2][1] else float("nan")
        line += f"  order {p:.2f}"
    print(line)
    ok &= max(abs(v[-1][1]), abs(v[-1][2])) < 1e-3 and (len(v) < 3 or p > 1.3)
print("PASS" if ok else "FAIL")
