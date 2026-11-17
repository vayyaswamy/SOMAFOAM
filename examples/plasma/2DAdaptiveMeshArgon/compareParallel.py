#!/usr/bin/env python3
"""Compare a serial and a parallel run of this case (paper case V8).

    python3 compareParallel.py SERIAL_DIR PARALLEL_DIR [--tol 1e-3]

The electrode voltage and total current written by the electrodeVoltageCurrent
function object (electrodes/0/electrode.dat) are compared at every step, and
the cell counts after each refinement are taken from log.somaFoam. Parallel
refinement keeps the 2:1 balance across processor boundaries and can select a
few more or fewer cells than the serial run, so the comparison uses a
tolerance on the current, relative to its largest magnitude, instead of
round-off. With the linear solvers converged to 1e-13 (Allrun), the runs
agreed to 4e-10 of the largest current before the first refinement and to
6e-5 after two changes of the mesh. The last line is PASS or FAIL.
"""
import argparse, os, re, sys


def read_dat(case):
    rows = []
    for line in open(os.path.join(case, "electrodes", "0", "electrode.dat")):
        if line.strip() and not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    return rows


def cells(case):
    log = os.path.join(case, "log.somaFoam")
    return [int(m.group(1)) for m in re.finditer(r"Mesh changed: (\d+) cells", open(log).read())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("serial")
    ap.add_argument("parallel")
    ap.add_argument("--tol", type=float, default=1e-3)
    a = ap.parse_args()
    s, p = read_dat(a.serial), read_dat(a.parallel)
    n = min(len(s), len(p))
    if n < 2 or abs(s[n - 1][0] - p[n - 1][0]) > 1e-3 * (s[1][0] - s[0][0]):
        print("the runs have different times\nFAIL")
        return 1
    imax = max(abs(r[4]) for r in s[:n])
    err = max(abs(s[i][4] - p[i][4]) for i in range(n)) / imax
    verr = max(abs(s[i][1] - p[i][1]) for i in range(n)) / max(abs(r[1]) for r in s[:n])
    print("steps compared            %d (t = %.4g s)" % (n - 1, s[n - 1][0]))
    print("cells after refinement    serial %s" % cells(a.serial))
    print("                          parallel %s" % cells(a.parallel))
    print("max |dV| / max |V|        %.3e" % verr)
    print("max |dI| / max |I|        %.3e (tolerance %.1e)" % (err, a.tol))
    ok = err <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
