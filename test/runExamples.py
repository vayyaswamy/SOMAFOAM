#!/usr/bin/env python3
"""Run the SOMAFOAM examples and collect their verdicts.

Every example directory under examples/ that has an Allrun script is a test. The
case is copied to a scratch directory, Allrun is run there, and the verdict is
the last line of its output that reads PASS or FAIL. A case without such a line
is reported as RAN (it ran to the end but has no reference yet) or ERROR (Allrun
returned a non-zero code or timed out).

    python3 test/runExamples.py                 # all examples
    python3 test/runExamples.py --list          # list the tests
    python3 test/runExamples.py flow icp        # examples whose path contains flow or icp
    python3 test/runExamples.py --short         # SOMAFOAM_TEST_SHORT=1 for Allrun
    python3 test/runExamples.py --jobs 4 --timeout 3600 --keep
    python3 test/runExamples.py --scratch /scratch/$USER/tests   # where the cases run

The source environment (etc/bashrc) must be loaded. The summary is printed and
written to test/results.json.
"""
import argparse, json, os, shutil, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")


def find_cases():
    cases = []
    for d, _, files in os.walk(EXAMPLES):
        if "Allrun" in files:
            cases.append(os.path.relpath(d, EXAMPLES))
    return sorted(cases)


def run_case(rel, args, scratch):
    src = os.path.join(EXAMPLES, rel)
    dst = os.path.join(scratch, rel)
    shutil.copytree(src, dst, symlinks=True)
    # cases may refer to shared scripts one level up (e.g. examples/flow/verify.py)
    parent_src, parent_dst = os.path.dirname(src), os.path.dirname(dst)
    for f in os.listdir(parent_src):
        p = os.path.join(parent_src, f)
        if os.path.isfile(p) and not os.path.exists(os.path.join(parent_dst, f)):
            shutil.copy2(p, parent_dst)
    env = dict(os.environ)
    if args.short:
        env["SOMAFOAM_TEST_SHORT"] = "1"
    t0 = time.time()
    try:
        r = subprocess.run(["sh", "./Allrun"], cwd=dst, env=env, capture_output=True,
                           text=True, timeout=args.timeout)
        out, code = r.stdout + r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out, code = (e.stdout or "") if isinstance(e.stdout, str) else "", "timeout"
    dt = time.time() - t0
    verdicts = [l.strip() for l in out.splitlines() if l.strip() in ("PASS", "FAIL")]
    if code != 0:
        status = "ERROR"
    elif verdicts:
        status = "FAIL" if "FAIL" in verdicts else "PASS"
    else:
        status = "RAN"
    with open(os.path.join(dst, "log.Allrun"), "w") as f:
        f.write(out)
    return {"case": rel, "status": status, "returncode": code, "seconds": round(dt, 1),
            "log": os.path.join(dst, "log.Allrun"), "tail": out.strip().splitlines()[-5:]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("patterns", nargs="*", help="run only cases whose path contains one of these")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--short", action="store_true", help="set SOMAFOAM_TEST_SHORT=1 for Allrun")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=7200.0, help="seconds per case")
    ap.add_argument("--keep", action="store_true", help="keep the scratch directory")
    ap.add_argument("--scratch", default=os.environ.get("SOMAFOAM_SCRATCH"),
                    help="parent directory of the test runs (default $SOMAFOAM_SCRATCH or /tmp)")
    args = ap.parse_args()

    cases = [c for c in find_cases() if not args.patterns or any(p in c for p in args.patterns)]
    if args.list:
        print("\n".join(cases)); return 0
    if not shutil.which("somaFoam"):
        print("somaFoam not found: load etc/bashrc first"); return 2

    if args.scratch:
        os.makedirs(args.scratch, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="somafoam_tests_", dir=args.scratch)
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        results = list(ex.map(lambda c: run_case(c, args, scratch), cases))

    width = max([len(r["case"]) for r in results] + [4])
    print(f"\n{'case':<{width}}  status  seconds")
    for r in results:
        print(f"{r['case']:<{width}}  {r['status']:<6}  {r['seconds']:>7}")
    n = {s: sum(r["status"] == s for r in results) for s in ("PASS", "FAIL", "RAN", "ERROR")}
    print(f"\n{n['PASS']} passed, {n['FAIL']} failed, {n['RAN']} ran without a reference, {n['ERROR']} errors")
    with open(os.path.join(ROOT, "test", "results.json"), "w") as f:
        json.dump({"short": args.short, "results": results}, f, indent=1)
    if args.keep:
        print("scratch kept in", scratch)
    else:
        shutil.rmtree(scratch, ignore_errors=True)
    return 1 if n["FAIL"] or n["ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
