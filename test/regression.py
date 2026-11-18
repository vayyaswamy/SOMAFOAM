"""Regression test of an example against stored reference fields.

    python3 test/regression.py store   <run directory> <reference directory> FIELD... [--regions R...]
    python3 test/regression.py compare <run directory> <reference directory> [--rtol 1e-6]

store writes the internal values of the given fields at the latest time of
the run (and of the same fields in the listed regions) to
reference/fields.json, at most 2000 values per field taken at an even
stride over the cells, so that the references of large meshes stay small.
compare reads the same values from the run and prints for each field the
largest difference relative to the largest value of the reference, and PASS
when all are below rtol. The
regression tests of the examples use it after a run of a fixed number of
steps (see the Allrun of plasma/1TorrArgonPlasma/defaultSchemes).
"""
import json
import os
import re
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


def sample(v):
    stride = max(1, -(-len(v) // 2000))
    return stride, v[::stride]


def store(run, ref, fields, regions):
    t = latest(run)
    out = {"time": os.path.basename(t), "fields": {}}
    for region in [""] + regions:
        for f in fields:
            path = os.path.join(t, region, f)
            if os.path.exists(path):
                v = values(path)
                stride, vs = sample(v)
                out["fields"][os.path.join(region, f)] = {"cells": len(v), "stride": stride, "values": vs}
                print("stored", os.path.join(region, f), len(vs), "of", len(v), "values from", out["time"])
    os.makedirs(ref, exist_ok=True)
    json.dump(out, open(os.path.join(ref, "fields.json"), "w"), indent=0)


def compare(run, ref, rtol):
    t = latest(run)
    data = json.load(open(os.path.join(ref, "fields.json")))
    ok = abs(float(os.path.basename(t)) / float(data["time"]) - 1) < 1e-9
    print(f"run time {os.path.basename(t)}, reference time {data['time']}{'' if ok else ' (differ)'}")
    for rel, d in sorted(data["fields"].items()):
        v = values(os.path.join(t, rel))
        b = d["values"]
        a = v[::d["stride"]] if len(v) == d["cells"] else []
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
