"""Rate tables of the argonGlobalModel case from the Maxwellian fits of the
argon global model of Lieberman and Lichtenberg (2005), with Te in volts:

    K_iz = 2.34e-14 Te^0.59 exp(-17.44/Te)                         m3/s
    K_ex = 2.48e-14 Te^0.33 exp(-12.78/Te)                         m3/s
    K_el = 2.336e-14 Te^1.609 exp(0.0618 ln^2 Te - 0.1171 ln^3 Te) m3/s

    python3 makeRates.py <constant directory>

writes reaction_iz, reaction_ex and reaction_el as (Te [K], ln k) with k in
m3/(kmol s), on 600 points from 0.2 to 30 eV, the format of the tabulated
electron-impact reactions. compare.py uses the same fits.
"""
import math
import os
import sys

AV, EV_K = 6.02214076e26, 11604.518


def Kiz(T):
    return 2.34e-14 * T**0.59 * math.exp(-17.44 / T)


def Kex(T):
    return 2.48e-14 * T**0.33 * math.exp(-12.78 / T)


def Kel(T):
    l = math.log(T)
    return 2.336e-14 * T**1.609 * math.exp(0.0618 * l * l - 0.1171 * l**3)


if __name__ == "__main__":
    d = sys.argv[1]
    head = "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      %s;\n}\n\n"
    for name, f in (("reaction_iz", Kiz), ("reaction_ex", Kex), ("reaction_el", Kel)):
        n = 600
        tv = [0.2 * (150.0)**(i / (n - 1)) for i in range(n)]     # 0.2 to 30 eV
        with open(os.path.join(d, name), "w") as fh:
            fh.write("// Lieberman and Lichtenberg (2005) Maxwellian fit, written by makeRates.py\n")
            fh.write("(\n")
            for t in tv:
                fh.write("(%.10g %.12g)\n" % (t * EV_K, math.log(max(f(t), 1e-60) * AV)))
            fh.write(")\n")
