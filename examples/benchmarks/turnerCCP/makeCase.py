"""Set up one case of the Turner et al. (2013) helium benchmark.

    python3 makeCase.py <case directory> <case 1-4> [periods]

The case directory is a copy of examples/plasma/HeliumIonDistribution, which
is case 4. This sets the gas density and pressure, the voltage, the uniform
initial density and the simulation time of table I of the benchmark (or the
given number of periods), and averages over the last 32 periods: the field
averages are reset at every write, and the fields are written every 32
periods of 250 steps.
"""
import re
import sys

KB, F = 1.380649e-23, 13.56e6
CASES = {  # N [m^-3], V [V], n0 [m^-3], periods
    1: (9.64e20, 450.0, 2.56e14, 1280),
    2: (3.21e21, 200.0, 5.12e14, 5120),
    3: (9.64e21, 150.0, 5.12e14, 5120),
    4: (3.21e22, 120.0, 3.84e14, 15360),
}
d, case = sys.argv[1], int(sys.argv[2])
N, V, n0, periods = CASES[case]
if len(sys.argv) > 3:
    periods = int(sys.argv[3])


def sub(path, pattern, repl, count=1):
    s = open(path).read()
    s2, n = re.subn(pattern, repl, s, count=count)
    if n == 0:
        raise SystemExit(f"{path}: no match for {pattern}")
    open(path, "w").write(s2)


sub(f"{d}/0/He", r"internalField\s+[^;]+;", f"internalField   uniform {N:.6e};")
sub(f"{d}/0/p", r"internalField\s+[^;]+;", f"internalField   uniform {N * KB * 300:.6e};")
sub(f"{d}/0/Phi", r"amplitude\s+uniform\s+[0-9.eE+-]+;", f"amplitude       uniform {V:g};")
sub(f"{d}/system/setExpressionFieldsDict", r"nWall\s+[0-9.eE+-]+;", f"nWall       {n0:.4e};")
sub(f"{d}/system/setExpressionFieldsDict", r"nCentre\s+[0-9.eE+-]+;", "nCentre     0;")
sub(f"{d}/system/controlDict", r"\nendTime\s+[^\n]*", f"\nendTime         {periods / F:.10e};    // {periods} periods")
sub(f"{d}/system/controlDict", r"\nwriteInterval\s+[^\n]*", "\nwriteInterval   8000;               // 32 periods of 250 steps")
# ion flux averaged as well, for the ion current at the electrodes
sub(f"{d}/system/controlDict", r"(\n\s*)N_Hep1(\s*\n\s*\{)", r"\1J_Hep1\2 mean on; prime2Mean off; base time; }\1N_Hep1\2")
print(f"case {case}: N = {N:.3e} m^-3, p = {N*KB*300:.4g} Pa, V = {V:g} V, n0 = {n0:.3e} m^-3, {periods} periods")
