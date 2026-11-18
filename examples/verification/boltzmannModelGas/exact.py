"""Model gases of the boltzmannModelGas case and their exact two-term solutions.

With elastic collisions only, the two-term equation of electronBoltzmann has
the zero-flux solution of Davydov,

    F0(eps) = C exp(-int_0^eps de / (kB T/e + (M/(6 m)) (E/N)^2/(e sigma(e)^2))),

eps in eV, E/N in V m2, normalised by int eps^1/2 F0 = 1. For a constant
cross section it is the Druyvesteyn distribution at T = 0, and for
sigma = sigma0 (eps0/eps)^1/2, a constant collision frequency, a Maxwellian of
kB Te/e = kB T/e + M (E/N)^2/(6 m sigma0^2 eps0). The mean energy is
int eps^3/2 F0 and the reduced mobility
mu N = -(gamma/3) int eps/sigma dF0/deps, gamma = (2 e/m)^1/2.

    python3 exact.py write <file>     writes the cross sections (LXCat format)
"""
import math
import sys

E, ME, KB = 1.602176634e-19, 9.1093837e-31, 1.380649e-23
MRATIO = 1.370558e-4            # m/M, helium
SIGMA0, T = 1e-19, 300.0
GASES = {"constantSigma": lambda e: SIGMA0, "constantNu": lambda e: SIGMA0 * math.sqrt(1.0 / e)}


def exact(gas, EN_Td, emax=None, n=400000):
    s = GASES[gas]
    EN = EN_Td * 1e-21
    kT = KB * T / E
    g = lambda e: 1.0 / (kT + EN * EN / (6 * MRATIO * e * s(e)**2))
    if emax is None:
        emax = 1.0
        while True:                      # until F0 is negligible
            de = emax / 2000
            integral = sum(g((i + 0.5) * de) * de for i in range(2000))
            if integral > 60:
                break
            emax *= 2
    de = emax / n
    eps, lnF, F = [], 0.0, []
    for i in range(n + 1):
        e = i * de
        if i:
            lnF -= g(e - 0.5 * de) * de
        eps.append(e)
        F.append(math.exp(lnF))
    # trapezoidal integrals
    def integ(vals):
        return de * (sum(vals) - 0.5 * (vals[0] + vals[-1]))
    norm = integ([math.sqrt(e) * f for e, f in zip(eps, F)])
    mean = integ([e**1.5 * f for e, f in zip(eps, F)]) / norm
    gamma = math.sqrt(2 * E / ME)
    dF = [(F[min(i + 1, n)] - F[max(i - 1, 0)]) / ((min(i + 1, n) - max(i - 1, 0)) * de) for i in range(n + 1)]
    muN = -gamma / 3 * integ([(e / s(max(e, 1e-12))) * d for e, d in zip(eps, dF)]) / norm
    return mean, muN


def write(path):
    with open(path, "w") as f:
        f.write("Model gases of examples/verification/boltzmannModelGas: elastic collisions only, mass ratio of helium\n\n")
        for gas, s in GASES.items():
            f.write("ELASTIC\n%s\n %.6e\nCOMMENT: %s\nCOLUMNS: Energy (eV) | Cross section (m2)\n-----------------------------\n"
                    % (gas, MRATIO, "constant cross section" if gas == "constantSigma" else "cross section 1e-19 (1 eV/eps)^1/2, constant collision frequency"))
            for i in range(2001):
                e = 1e-4 * (1e7)**(i / 2000)        # 1e-4 to 1000 eV
                f.write(" %.6e\t%.6e\n" % (e, s(e)))
            f.write("-----------------------------\n\n")


if __name__ == "__main__":
    if sys.argv[1] == "write":
        write(sys.argv[2])
    else:
        for gas in GASES:
            for en in (0.5, 1, 2, 4):
                print(gas, en, "Td: mean energy %.6f eV, mu N %.6e" % exact(gas, en))
