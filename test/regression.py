"""Regression test of an example against stored reference fields.

    python3 test/regression.py store   <run directory> <reference directory> FIELD... [--regions R...]
    python3 test/regression.py compare <run directory> <reference directory> [--rtol 1e-6]

store copies the internal values of the given fields at the latest time of
the run (and of the same fields in the listed regions) to the reference
directory. compare reads every field in the reference directory and its
regions, and prints for each the largest difference relative to the
largest value of the reference, and PASS when all are below rtol. The
regression tests of the examples use it after a run of a fixed number of
steps (see the Allrun of plasma/1TorrArgonPlasma/defaultSchemes).
"""
import os
import re
import shutil
import sys


def latest(run):
    times = [d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d) and float(d) > 0]
    if not times:
        sys.exit(f"{run}: no time directory after 0")
    return os.path.join(run, max(times, key=float))


def values(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    if m:
        n = int(m.group(2))
        body = s[m.end():]
        if m.group(1) == "scalar":
            return [float(x) for x in body.split(")")[0].split()][:n]
        return [float(x) for v in re.findall(r"\(([^()]*)\)", body)[:n] for x in v.split()]
    m = re.search(r"internalField\s+uniform\s+(\([^)]*\)|\S+);", s)
    return [float(x) for x in m.group(1).strip("()").split()]


def store(run, ref, fields, regions):
    t = latest(run)
    for region in [""] + regions:
        os.makedirs(os.path.join(ref, region), exist_ok=True)
        for f in fields:
            src = os.path.join(t, region, f)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(ref, region, f))
                print("stored", os.path.join(region, f), "from", os.path.basename(t))
    open(os.path.join(ref, "TIME"), "w").write(os.path.basename(t) + "\n")


def compare(run, ref, rtol):
    t = latest(run)
    want = open(os.path.join(ref, "TIME")).read().split()[0]
    ok = abs(float(os.path.basename(t)) / float(want) - 1) < 1e-9
    print(f"run time {os.path.basename(t)}, reference time {want}{'' if ok else ' (differ)'}")
    for root, _, files in os.walk(ref):
        for f in sorted(files):
            if f == "TIME":
                continue
            rel = os.path.relpath(os.path.join(root, f), ref)
            a, b = values(os.path.join(t, rel)), values(os.path.join(ref, rel))
            scale = max(abs(x) for x in b) or 1.0
            err = max(abs(x - y) for x, y in zip(a, b)) / scale if len(a) == len(b) else float("inf")
            good = err < rtol
            ok &= good
            print(f"  {rel:24s} max difference {err:.2e} of max |ref| {scale:.4g}  {'ok' if good else 'DIFFERS'}")
    print("PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    mode, run, ref = sys.argv[1:4]
    rest = sys.argv[4:]
    if mode == "store":
        regions = rest[rest.index("--regions") + 1:] if "--regions" in rest else []
        fields = rest[:rest.index("--regions")] if "--regions" in rest else rest
        store(run, ref, fields, regions)
    else:
        rtol = float(rest[rest.index("--rtol") + 1]) if "--rtol" in rest else 1e-6
        sys.exit(0 if compare(run, ref, rtol) else 1)
