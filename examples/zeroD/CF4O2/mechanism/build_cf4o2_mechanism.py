#!/usr/bin/env python3
"""CF4/O2 plasma mechanism for SOMAFOAM.

Writes the chemistry of a SOMAFOAM case: constant/reactions,
constant/thermoProperties, the electron transport tables mu_electron and
D_electron, the Te tables of the electron-impact processes that come from the
Boltzmann equation, and MECHANISM.md, a table of every reaction with its rate,
energy loss and source.

Usage: build_cf4o2_mechanism.py <boltzmannDir> <caseDir> [p_Pa] [Tg_K]
  boltzmannDir  electronBoltzmann output (constant/boltzmann of the CF4/O2
                case, 80 % CF4 and 20 % O2, LXCat Hayashi CF4 and Phelps O2)
                with the cross-section file in ../crossSections
  caseDir       SOMAFOAM case to write into (constant/ must exist)
  p_Pa, Tg_K    gas pressure and temperature, 133.322 Pa (1 Torr) and 300 K
                by default

The mechanism
-------------
Base set: Kimura and Noto, J. Appl. Phys. 100, 063303 (2006), Table I, the
reactions k1-k64 with their species, rate coefficients and products, written
as published (m3/s, Te in eV), with these changes.

1. Association reactions at the operating pressure. Kimura's k6, k7 and k11
   (CF3 + F, CF2 + F, CF + F), which are linear in p in mTorr, are replaced
   by the fall-off rates of Plumb and Ryan, Plasma Chem. Plasma Process. 6,
   205 (1986), Table II (k0, kinf, Fc for CF4 as the third body):
   k = k0 M / (1 + k0 M/kinf) Fc^(1/(1 + log10(k0 M/kinf)^2)).
   Plumb and Ryan's starred reactions that Kimura does not have are added:
   COF + F -> COF2 and F + CO -> COF (fall-off, Table II) and
   CF2 + O -> CO + 2F (Table I). The three-body reactions F + F + M (k24) and
   F + O + M (k53) of Kimura and O + O2 + M -> O3 + M (Atkinson et al.,
   Atmos. Chem. Phys. 4, 1461 (2004), M = O2) are multiplied by the total
   density M = p/(kB Tg). Pressure and Tg are fixed in the cases of this
   work, so these are the exact rates there.
2. Species that Kimura produces but does not carry are added: O3 (k47),
   O2- (k46) and C (k61). O3, O2- and O2(a1Dg) get the heavy-particle
   reactions of Toneli et al., J. Phys. D 48, 325202 (2015), Tables 2-4,
   for the species of this set; reactions that Kimura already has keep
   Kimura's rate. O(1D) is lumped into O, as in Kimura. O3- is not carried,
   so Toneli's reactions that form it (99, 108) are left out.
   O- + O2 -> O2- + O and O- + O2 -> O + O2 + e (Toneli 112, 113) are left
   out: they need 1.0 and 1.46 eV, and their constant rates hold for fast
   ions, not for ions at the gas temperature.
3. Ion-ion recombination. Kimura's X+ + F- -> X + F (k18) and
   X+ + O- -> X + O (k33, k32 for O+) are applied to every positive ion of the
   set. For O2-, Toneli 105-107 are used with O2+ and O+, and the sum of
   Toneli 105 and 106 is used for the other positive ions, for which no
   O2- rate was found.
4. Electron impact on the feed gas. Kimura's fits assume a Maxwellian
   distribution, which in CF4 has far more fast electrons than the
   distribution of the Boltzmann equation; at the same Te his CF4 ionization
   rate was 10-100 times that of the cross sections between 2 and 4 eV, and
   the net ionization of the feed gas balanced at about 61 Td instead of
   137 Td. The rates of CF4 and O2 are therefore those of the LXCat cross
   sections (Hayashi CF4, Phelps O2) over the distributions of the Boltzmann
   equation for the feed gas (SOMAFOAM electronBoltzmann, two-term), tabulated
   in Te, with these channels.
   CF4 neutral dissociation (Hayashi 12.85 eV, qd - qi) split into CF3 + F,
   CF2 + 2F and CF + 3F (k1-k3) with the fragment fractions of Christophorou
   et al., J. Phys. Chem. Ref. Data 25, 1341 (1996), Table 21 (all CF3 below
   16 eV). CF4 total ionization split into CF3+ + F, CF2+ + 2F, CF+ + 3F
   (k8-k10), C+ + 4F and F+ + CF3 with the partial ionization cross sections
   of their Table 15; the CF3+ entries printed there as 2.66 (55 eV) and 2.95
   (70 eV) are read as 3.66 and 3.95, the averages of Tables 13 and 14 that
   the table is made of. CF4 total attachment (F- and CF3-) as F- + CF3
   (k15). O2 ionization to O2+ (k30), dissociative attachment (k31), the
   6.0 and 8.4 eV losses as dissociation into 2 O (k35, k36), O2(a1Dg)
   excitation (k41), and its de-excitation (k44) by detailed balance,
   sigma(e) = (3/2) (e + u)/e sigma_a1(e + u). Kimura's k37 (O- + O+) and k38
   (O+ + O) are left out, because the Phelps ionization cross section is the
   total one and its channels are not given. The three-body attachment of the
   Phelps set is left out; with O2 as third body it is more than three orders
   below the CF4 attachment. Every channel carries the threshold of its cross
   section as energy loss, so that the electron energy balance is the one of
   the Boltzmann solution; the channels of each process add up to the rate of
   electronBoltzmann, which the script checks.
5. The remaining processes of the two sets (elastic, vibrational, CF4 7.86
   eV, O2 b1, 4.5 and 9.97 eV) are energy-loss-only collisions with the
   Boltzmann rate. Electron transport comes from the same solution.
6. Electron impact on the radicals and products keeps Kimura's fits, with the
   activation energy of the fit as energy loss. Attachment and recombination
   carry no loss.
"""
import math
import pathlib
import re
import sys

import numpy as np

AV = 6.02214076e26          # 1/kmol, as in SOMAFOAM
KB = 1.380649e-23
QE = 1.602176634e-19
ME = 9.1093837e-31
EV_K = QE / KB              # 11604.518 K per eV
WE = 5.48579909e-4          # electron, kg/kmol

bdir = pathlib.Path(sys.argv[1])
case = pathlib.Path(sys.argv[2])
P = float(sys.argv[3]) if len(sys.argv) > 3 else 133.322
TG = float(sys.argv[4]) if len(sys.argv) > 4 else 300.0
M = P / (KB * TG)            # total density, m^-3
const = case / "constant"

# --------------------------------------------------------------------------
# Species: name, composition, charge. Charged species first (activeSpecies).
AT = {"C": 12.011, "F": 18.998403, "O": 15.999}
SPECIES = [
    ("electron", {}, -1),
    ("CF3p1", {"C": 1, "F": 3}, 1), ("CF2p1", {"C": 1, "F": 2}, 1),
    ("CFp1", {"C": 1, "F": 1}, 1), ("Fp1", {"F": 1}, 1),
    ("O2p1", {"O": 2}, 1), ("Op1", {"O": 1}, 1),
    ("COp1", {"C": 1, "O": 1}, 1), ("CO2p1", {"C": 1, "O": 2}, 1),
    ("Cp1", {"C": 1}, 1),
    ("Fn1", {"F": 1}, -1), ("On1", {"O": 1}, -1), ("O2n1", {"O": 2}, -1),
    ("CF3", {"C": 1, "F": 3}, 0), ("CF2", {"C": 1, "F": 2}, 0),
    ("CF", {"C": 1, "F": 1}, 0), ("F", {"F": 1}, 0), ("F2", {"F": 2}, 0),
    ("O2", {"O": 2}, 0), ("O2a", {"O": 2}, 0), ("O", {"O": 1}, 0),
    ("O3", {"O": 3}, 0), ("CO", {"C": 1, "O": 1}, 0),
    ("CO2", {"C": 1, "O": 2}, 0), ("FO", {"F": 1, "O": 1}, 0),
    ("COF", {"C": 1, "O": 1, "F": 1}, 0),
    ("COF2", {"C": 1, "O": 1, "F": 2}, 0), ("C", {"C": 1}, 0),
    ("CF4", {"C": 1, "F": 4}, 0),
]
NAMES = [s[0] for s in SPECIES]
COMP = {s[0]: s[1] for s in SPECIES}
CHARGE = {s[0]: s[2] for s in SPECIES}
N_CHARGED = sum(1 for s in SPECIES if s[2] != 0)
assert all(s[2] != 0 for s in SPECIES[:N_CHARGED])


def molar_mass(name):
    if name == "electron":
        return WE
    return sum(AT[e] * n for e, n in COMP[name].items()) - CHARGE[name] * WE


POS = [s for s in NAMES if CHARGE[s] > 0]
NEUTRAL_OF = {"CF3p1": "CF3", "CF2p1": "CF2", "CFp1": "CF", "Fp1": "F",
              "O2p1": "O2", "Op1": "O", "COp1": "CO", "CO2p1": "CO2",
              "Cp1": "C"}

# --------------------------------------------------------------------------
# Boltzmann results of the feed gas
proc = []   # (index, type, target line, parameter)
for line in (bdir.parent.parent / "log.electronBoltzmann").read_text().splitlines():
    m = re.match(r"^\s+(\d+)\s+(ATTACHMENT|ELASTIC|EFFECTIVE|EXCITATION|"
                 r"IONIZATION)\s+(.*?)\s+parameter\s+(\S+),", line)
    if m:
        proc.append((int(m.group(1)), m.group(2), m.group(3),
                     float(m.group(4))))


def proc_index(kind, target):
    hits = [p[0] for p in proc if p[1] == kind and p[2] == target]
    assert len(hits) == 1, (kind, target, hits)
    return hits[0]


def read_table(path):
    rows = re.findall(r"^\s*\(?\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s*\)?\s*$",
                      path.read_text(), re.M)
    return np.array([[float(a), float(b)] for a, b in rows])


# Te of each E/N and the distributions, for the O2 elastic rate
tr = np.loadtxt(bdir / "transport.dat", comments="#")
EN, TE = tr[:, 0], tr[:, 2]
bad = np.where(np.diff(TE) <= 0)[0]
assert bad.size == 0, ("Te of the Boltzmann solution not increasing with E/N at "
                       f"{EN[bad + 1]} Td; refine nEnergyCells")


def lxcat_blocks(path):
    """Cross sections of an LXCat file: list of (type, target, param, E, s)."""
    out = []
    lines = path.read_text().replace("\r", "").splitlines()
    i = 0
    while i < len(lines):
        k = lines[i].strip()
        if k in ("ELASTIC", "EFFECTIVE", "EXCITATION", "IONIZATION",
                 "ATTACHMENT"):
            target = lines[i + 1].strip()
            par = 0.0 if k == "ATTACHMENT" else float(lines[i + 2].split()[0])
            j = i + 1
            while not lines[j].startswith("-----"):
                j += 1
            j += 1
            e, s = [], []
            while not lines[j].startswith("-----"):
                a, b = lines[j].split()[:2]
                e.append(float(a))
                s.append(float(b))
                j += 1
            out.append((k, target, par, np.array(e), np.array(s)))
            i = j
        i += 1
    return out


def sigma(block, eps):
    """As electronBoltzmann: linear, zero below a threshold process, else
    constant beyond the ends."""
    k, _, _, e, s = block
    out = np.interp(eps, e, s)
    if k not in ("ELASTIC", "EFFECTIVE"):
        out = np.where(eps <= e[0], 0.0, out)
    return out


xsfile = next((bdir.parent / "crossSections").glob("*.txt"))
blocks = [b for b in lxcat_blocks(xsfile)
          if b[1].split()[0] in ("CF4", "O2")]


DIST = []
_files = list(bdir.glob("distribution_*.dat"))
for en in EN:
    f = min(_files, key=lambda p: abs(float(p.stem.split("_", 1)[1]) - en))
    d = np.loadtxt(f, comments="#")
    DIST.append((d[:, 0], d[:, 1], np.gradient(d[:, 0])))


def rate_from_F0(sig_fn):
    """k(Te) = sqrt(2e/m) int eps sigma F0 d eps over the distributions."""
    return np.array([math.sqrt(2 * QE / ME) * np.sum(eps * sig_fn(eps) * F0 * de)
                     for eps, F0, de in DIST])


def block(kind, target):
    hits = [b for b in blocks if b[0] == kind and b[1] == target]
    assert len(hits) == 1, (kind, target)
    return hits[0]


def check_total(kind, target, k):
    """The channels of a process must add up to the rate of the solver."""
    ref = tr[:, 7 + proc_index(kind, target)]
    s = ref > 1e-22
    dev = np.max(np.abs(k[s] / ref[s] - 1)) if s.any() else 0.0
    print(f"{kind} {target}: channels vs electronBoltzmann, max relative "
          f"difference {dev:.2e}")
    assert dev < 1e-3, (kind, target, dev)


o2 = [b for b in blocks if b[1].split()[0] == "O2"]
o2_eff = [b for b in o2 if b[0] == "EFFECTIVE"][0]
o2_inel = [b for b in o2 if b[0] in ("EXCITATION", "IONIZATION")]
k_o2_el = rate_from_F0(lambda e: np.maximum(
    sigma(o2_eff, e) - sum(sigma(b, e) for b in o2_inel), 0.0))
# Check of the integration: the momentum-transfer rate of the solver
k_o2_eff = rate_from_F0(lambda e: sigma(o2_eff, e))
i_eff = proc_index("EFFECTIVE", "O2")
dev = np.max(np.abs(k_o2_eff / tr[:, 7 + i_eff] - 1))
print(f"O2 effective rate from the distributions vs electronBoltzmann: "
      f"max relative difference {dev:.2e}")
assert dev < 2e-2, dev

# --------------------------------------------------------------------------
# Reactions
R = []      # dicts: type, lhs, rhs, args, dE (J), note, src


def side(s):
    out = []
    for t in s.split("+"):
        t = t.strip()
        m = re.match(r"^(\d+)\s*(\S+)$", t)
        n, name = (int(m.group(1)), m.group(2)) if m else (1, t)
        assert name in NAMES, name
        out += [name] * n
    return out


def e_arr(label, eq, A, b, E, dE, src):
    """Electron impact, k = A Te^b exp(-E/Te), m3/s, Te in eV."""
    lhs, rhs = (side(x) for x in eq.split("->"))
    assert lhs[0] == "electron"
    R.append(dict(type="electronImpactInelasticArrheniusReaction", lhs=lhs,
                  rhs=rhs, A=A * AV * EV_K ** (-b), b=b, Ta=E * EV_K,
                  dE=dE * QE, label=label, src=src,
                  rate=f"{A:.3g} Te^{b:g} exp(-{E:g}/Te)" if b else
                  f"{A:.3g} exp(-{E:g}/Te)" if E else f"{A:.3g}",
                  dEeV=dE))


def e_tab(label, eq, table, dE, src, elastic=False):
    lhs, rhs = (side(x) for x in eq.split("->"))
    R.append(dict(type=("electronImpactElasticTabularReaction" if elastic
                        else "electronImpactInelasticTabularReaction"),
                  lhs=lhs, rhs=rhs, table=table, dE=dE * QE, label=label,
                  src=src, rate="Boltzmann table", dEeV=dE))


def heavy(label, eq, k, src, b=0.0, Ta=0.0, order=None, rate=None):
    """k (Tg/300)^b exp(-Ta/Tg) in SI (1/s, m3/s)."""
    lhs, rhs = (side(x) for x in eq.split("->"))
    n = len(lhs) if order is None else order
    A = k * 300.0 ** (-b) * AV ** (n - 1)
    R.append(dict(type="irreversibleArrheniusReaction", lhs=lhs, rhs=rhs,
                  A=A, b=b, Ta=Ta, dE=0.0, label=label, src=src,
                  rate=rate or (f"{k:.3g}" + (f" (Tg/300)^{b:g}" if b else "")
                                + (f" exp(-{Ta:g}/Tg)" if Ta else "")),
                  dEeV=None))


def falloff(k0, kinf, Fc, Mc):
    """Plumb and Ryan Table II; cm units, Mc in cm^-3; returns cm3/s."""
    x = k0 * Mc / kinf
    return k0 * Mc / (1 + x) * Fc ** (1.0 / (1 + math.log10(x) ** 2))


# Check against Plumb and Ryan Table I (rates at 0.5 Torr, M = 1.6e16 cm-3)
PR = {6: (7.7e-27, 2e-11, 0.63, 1.3e-11), 8: (3.0e-29, 1.3e-11, 0.73, 4.2e-13),
      9: (3.2e-31, 1e-11, 0.72, 5e-15), 31: (6.5e-29, 1.4e-11, 0.68, 8e-13),
      38: (8.1e-32, 9.4e-11, 0.73, 1.3e-15)}
for n, (k0, ki, fc, tab) in PR.items():
    k = falloff(k0, ki, fc, 1.6e16)
    print(f"Plumb and Ryan reaction {n}: 0.5 Torr {k:.3g} cm3/s, Table I {tab:.3g}")
    assert abs(k / tab - 1) < 0.08, (n, k, tab)
Mc = M * 1e-6


def pr(n):
    k0, ki, fc, _ = PR[n]
    return falloff(k0, ki, fc, Mc) * 1e-6     # m3/s


KN = "Kimura and Noto (2006), Table I"
TH_FIT = "activation energy of the fit"
tables = {}      # Te tables written to constant/

# Electron impact on CF4 and O2: rates of the LXCat cross sections over the
# Boltzmann distributions (tabulated in Te), split into Kimura's channels
def e_xs(label, eq, sig_fn, dE, src):
    k = rate_from_F0(sig_fn)
    name = f"rate_{label}"
    tables[name] = np.column_stack([TE, np.log(np.maximum(k * AV, 1e-50 * AV))])
    e_tab(label, eq, name, dE, src)
    R[-1]["rate"] = "cross section over the EEDF"
    return k


def branch(e_tab_, cols):
    """Fractions of the columns of a table of partial cross sections; all of
    the first column below the first energy, constant above the last."""
    e = np.array([r[0] for r in e_tab_])
    v = np.array([r[1:] for r in e_tab_])
    tot = v.sum(axis=1)

    def frac(j):
        def f(eps):
            num = np.interp(eps, e, v[:, j])
            den = np.interp(eps, e, tot)
            out = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
            return np.where(eps < e[0], 1.0 if j == 0 else 0.0, out)
        return f
    return [frac(j) for j in range(len(cols))]


# Christophorou et al. (1996), Table 15: partial ionization cross sections of
# CF4 for CF3+, CF2+, CF+, C+, F+ (1e-20 m2). The CF3+ entries printed as 2.66
# (55 eV) and 2.95 (70 eV) are 3.66 and 3.95, the averages of Tables 13 and 14
# that Table 15 is made of (3.81 and 3.51, 4.13 and 3.77).
T15 = [(16, .02, 0, 0, 0, 0), (18, .02, 0, 0, 0, 0), (20, .31, 0, 0, 0, 0),
       (25, 1.04, .03, 0, 0, 0), (30, 1.89, .13, 0, 0, 0), (35, 2.41, .23, .03, 0, 0),
       (40, 2.90, .23, .12, .03, .02), (45, 3.23, .27, .21, .10, .08),
       (50, 3.48, .30, .25, .15, .14), (55, 3.66, .31, .28, .17, .17),
       (60, 3.81, .33, .31, .21, .23), (65, 3.88, .34, .34, .24, .29),
       (70, 3.95, .34, .36, .26, .32), (75, 3.99, .35, .36, .27, .33),
       (80, 4.03, .35, .36, .28, .34), (85, 4.06, .36, .36, .28, .35),
       (90, 4.08, .36, .36, .29, .36), (95, 4.09, .35, .36, .29, .39),
       (100, 4.10, .35, .36, .30, .39), (110, 4.13, .35, .36, .30, .39),
       (120, 4.12, .35, .35, .29, .39), (130, 4.08, .34, .33, .28, .39),
       (140, 4.04, .34, .32, .28, .39), (150, 4.00, .34, .30, .28, .39),
       (160, 3.95, .33, .29, .27, .38), (170, 3.90, .33, .28, .27, .38),
       (180, 3.86, .32, .27, .25, .35), (190, 3.83, .32, .26, .25, .35),
       (200, 3.79, .31, .26, .25, .34), (250, 3.62, .30, .24, .23, .33),
       (300, 3.48, .28, .21, .20, .26), (350, 3.36, .27, .19, .19, .24),
       (400, 3.25, .26, .19, .19, .25), (450, 3.16, .25, .18, .18, .23),
       (500, 3.08, .25, .17, .17, .23)]
assert abs(3.66 - (3.81 + 3.51) / 2) < 6e-3 and abs(3.95 - (4.13 + 3.77) / 2) < 6e-3
# Christophorou et al. (1996), Table 21: neutral CF3, CF2, CF fragments of CF4
# (Nakano and Sugai as revised by Sugai et al.), 1e-20 m2
T21 = [(16, .001, 0, 0), (18, .002, .0002, 0), (30, .015, .0017, .0012),
       (40, .020, .0031, .0025), (50, .032, .0042, .0042), (60, .044, .0057, .0069),
       (70, .052, .0066, .0105), (80, .065, .0081, .0128), (90, .067, .0095, .0148),
       (100, .074, .0108, .0152), (110, .089, .0114, .0152), (120, .093, .0118, .0152),
       (130, .093, .0120, .0152), (140, .091, .0119, .0150), (150, .089, .0116, .0148),
       (160, .087, .0113, .0147), (170, .085, .0110, .0143), (180, .084, .0106, .0141),
       (190, .081, .0102, .0132), (200, .077, .0101, .0134), (210, .076, .0097, .0127),
       (220, .074, .0096, .0123), (230, .073, .0091, .0122), (240, .073, .0086, .0119),
       (250, .071, .0084, .0114), (260, .070, .0082, .0108), (270, .068, .0080, .0105),
       (280, .066, .0076, .0105), (290, .066, .0076, .0100), (300, .065, .0072, .0097)]
# Check of the transcription against the sums of Table 22
for e_, s22 in ((30, .018), (50, .040), (100, .100), (150, .115), (200, .101), (300, .082)):
    r_ = [r for r in T21 if r[0] == e_][0]
    assert abs(sum(r_[1:]) - s22) < 1.5e-3, (e_, sum(r_[1:]), s22)

HY = "LXCat Hayashi CF4 "
CH = "Christophorou et al. (1996)"
s_dis = block("EXCITATION", "CF4 -> CF4*(12.85eV)")
s_ion = block("IONIZATION", "CF4 -> CF4^+")
s_att = block("ATTACHMENT", "CF4 -> CF4^-")
u_dis, u_ion = s_dis[2], s_ion[2]
fd = branch(T21, ("CF3", "CF2", "CF"))
k_sum = 0
for j, (lab, eq) in enumerate((("k1", "electron + CF4 -> electron + CF3 + F"),
                               ("k2", "electron + CF4 -> electron + CF2 + 2F"),
                               ("k3", "electron + CF4 -> electron + CF + 3F"))):
    k_sum = k_sum + e_xs(lab, eq, lambda e, j=j: sigma(s_dis, e) * fd[j](e), u_dis,
                         HY + "neutral dissociation (12.85 eV) times the CF" + "3 2 1".split()[j]
                         + " fraction of " + CH + ", Table 21 (replaces Kimura " + lab + ")")
check_total("EXCITATION", "CF4 -> CF4*(12.85eV)", k_sum)
fi = branch(T15, ("CF3+", "CF2+", "CF+", "C+", "F+"))
k_sum = 0
for j, (lab, eq) in enumerate((("k8", "electron + CF4 -> 2electron + CF3p1 + F"),
                               ("k9", "electron + CF4 -> 2electron + CF2p1 + 2F"),
                               ("k10", "electron + CF4 -> 2electron + CFp1 + 3F"),
                               ("CXi", "electron + CF4 -> 2electron + Cp1 + 4F"),
                               ("FXi", "electron + CF4 -> 2electron + Fp1 + CF3"))):
    k_sum = k_sum + e_xs(lab, eq, lambda e, j=j: sigma(s_ion, e) * fi[j](e), u_ion,
                         HY + "total ionization times the " + ("CF3+", "CF2+", "CF+", "C+", "F+")[j]
                         + " fraction of " + CH + ", Table 15"
                         + (f" (replaces Kimura {lab})" if lab.startswith("k") else ""))
check_total("IONIZATION", "CF4 -> CF4^+", k_sum)
k = e_xs("k15", "electron + CF4 -> Fn1 + CF3", lambda e: sigma(s_att, e), 0.0,
         HY + "total attachment (F- and CF3-) as F- + CF3 (replaces Kimura k15)")
check_total("ATTACHMENT", "CF4 -> CF4^-", k)
PH = "LXCat Phelps O2 "
for lab, eq, tgt, kind in (
        ("k30", "electron + O2 -> 2electron + O2p1", "O2 -> O2^+", "IONIZATION"),
        ("k31", "electron + O2 -> O + On1", "O2 -> O^-+O", "ATTACHMENT"),
        ("k35", "electron + O2 -> electron + 2O", "O2 -> O2(6.0eV)", "EXCITATION"),
        ("k36", "electron + O2 -> electron + 2O", "O2 -> O2(8.4eV)", "EXCITATION"),
        ("k41", "electron + O2 -> electron + O2a", "O2 -> O2(a1)", "EXCITATION")):
    b_ = block(kind, tgt)
    k = e_xs(lab, eq, lambda e, b_=b_: sigma(b_, e), 0.0 if kind == "ATTACHMENT" else b_[2],
             PH + tgt + f" (replaces Kimura {lab})")
    check_total(kind, tgt, k)
# De-excitation of O2(a1Dg) by detailed balance with the a1 cross section,
# sigma_sup(e) = (g_X/g_a) (e + u)/e sigma(e + u), g_X = 3, g_a = 2
b_a = block("EXCITATION", "O2 -> O2(a1)")
u_a = b_a[2]
e_xs("k44", "electron + O2a -> electron + O2",
     lambda e: 1.5 * (e + u_a) / np.maximum(e, 1e-12) * sigma(b_a, e + u_a), -u_a,
     PH + "O2(a1) by detailed balance (replaces Kimura k44)")
# Three-body attachment e + O2 + M -> O2- + M of the Phelps set (normalized to
# 1 cm-3) is left out: with O2 as third body at the operating pressure its
# rate is far below that of the dissociative attachment of CF4
b3 = block("ATTACHMENT", "O2")
k3b = rate_from_F0(lambda e: sigma(b3, e)) * 0.2 * M * 1e-6
k_cf4att = rate_from_F0(lambda e: sigma(s_att, e))
j100 = np.argmin(np.abs(EN - 100))
print(f"three-body O2 attachment at 100 Td {k3b[j100]:.2e} m3/s, CF4 attachment "
      f"{k_cf4att[j100]:.2e} m3/s")
assert k3b[j100] < 1e-3 * k_cf4att[j100]
# Electron impact on the radicals and products (Kimura, Maxwellian fits)
e_arr("k4", "electron + CF3 -> electron + CF2 + F", 6.48e-14, -0.96, 11.25, 11.25, KN + "; dE " + TH_FIT)
e_arr("k5", "electron + CF2 -> electron + CF + F", 1.16e-14, -0.38, 14.35, 14.35, KN + "; dE " + TH_FIT)
e_arr("k12", "electron + CF3 -> 2electron + CF3p1", 8e-15, 0, 12.2, 12.2, KN + "; dE " + TH_FIT)
e_arr("k13", "electron + CF2 -> 2electron + CF2p1", 2.5e-14, 0, 12.2, 12.2, KN + "; dE " + TH_FIT)
e_arr("k14", "electron + CF -> 2electron + CFp1", 2.5e-14, 0, 15.3, 15.3, KN + "; dE " + TH_FIT)
e_arr("k16a", "electron + CF2p1 -> CF2", 6e-14, 0, 0, 0.0, KN)
e_arr("k16b", "electron + CF3p1 -> CF3", 6e-14, 0, 0, 0.0, KN)
e_arr("k19", "electron + CF3 -> 2electron + CF2p1 + F", 2.20e-14, 0, 18.6, 18.6, KN + "; dE " + TH_FIT)
e_arr("k20", "electron + CF2 -> 2electron + CFp1 + F", 2.33e-14, 0, 17.0, 17.0, KN + "; dE " + TH_FIT)
e_arr("k21", "electron + CF3 -> CF2 + Fn1", 1.6e-16, -0.4329, 0.18, 0.0, KN)
e_arr("k22", "electron + CF2 -> CF + Fn1", 1.6e-16, -0.4329, 1.93, 0.0, KN)
e_arr("k23", "electron + F -> 2electron + Fp1", 1.3e-14, 0, 16.5, 16.5, KN + "; dE " + TH_FIT)
e_arr("k25", "electron + F2 -> Fn1 + F", 4.5e-15, -1.35, 0.15, 0.0, KN)
e_arr("k26", "electron + F2 -> electron + 2F", 1.18e-14, 0, 5.77, 5.77, KN + "; dE " + TH_FIT)
e_arr("k29", "electron + O -> 2electron + Op1", 9.0e-15, 0.7, 13.6, 13.6, KN + "; dE " + TH_FIT)
e_arr("k34", "electron + On1 -> 2electron + O", 2.0e-13, 0, 5.5, 5.5, KN + "; dE " + TH_FIT)
e_arr("k39", "electron + O2p1 -> 2O", 5.2e-15, -1, 0, 0.0, KN)
e_arr("k42", "electron + O2a -> 2electron + O2p1", 9.0e-16, 2, 11.6, 11.6, KN + "; dE " + TH_FIT)
e_arr("k43", "electron + O2a -> On1 + O", 2.28e-16, 0, 2.29, 0.0, KN)
e_arr("k45", "electron + O2a -> electron + 2O", 4.2e-15, 2, 4.6, 4.6, KN + "; dE " + TH_FIT)
# Electron impact on COF2, CO2, CO
e_arr("k59", "electron + COF2 -> electron + COF + F", 6.48e-14, -0.96, 11.25, 11.25, KN + " (k59 = k4); dE " + TH_FIT)
e_arr("k60", "electron + CO2 -> electron + CO + O", 8.5e-16, 0, 0, 0.0, KN + " (constant)")
e_arr("k61", "electron + CO -> electron + C + O", 1.87e-14, 0, 13.89, 13.89, KN + "; dE " + TH_FIT)
e_arr("k62", "electron + CO2 -> 2electron + CO2p1", 5.03e-14, 0, 16.26, 16.26, KN + "; dE " + TH_FIT)
e_arr("k63", "electron + CO -> 2electron + COp1", 4.93e-14, 0, 16.76, 16.76, KN + "; dE " + TH_FIT)

# Energy-loss-only collisions with the feed gas (Boltzmann tables)
def bz(kind, target):
    i = proc_index(kind, target)
    t = read_table(bdir / f"rate_{i}")
    name = f"rate_{i}"
    tables[name] = t
    return name


BZ = "Boltzmann equation, LXCat "
e_tab("elCF4", "electron + CF4 -> electron + CF4", bz("ELASTIC", "CF4"), 0.0, BZ + "Hayashi CF4", elastic=True)
tables["rate_O2_elastic"] = np.column_stack([TE, np.log(np.maximum(k_o2_el * AV, 1e-50 * AV))])
e_tab("elO2", "electron + O2 -> electron + O2", "rate_O2_elastic", 0.0, BZ + "Phelps O2, effective minus inelastic", elastic=True)
for tgt, dE in (("CF4 -> CF4(V24)", 0.052), ("CF4 -> CF4(V13)", 0.109),
                ("CF4 -> CF4*(7.86eV)", 7.86)):
    e_tab(tgt.split("(")[1].rstrip(")") if "V" in tgt else "CF4 7.86",
          "electron + CF4 -> electron + CF4", bz("EXCITATION", tgt), dE, BZ + "Hayashi CF4, " + tgt)
for tgt in ("O2 -> O2(rot)", "O2 -> O2(v1)", "O2 -> O2(v1res)", "O2 -> O2(v2)",
            "O2 -> O2(v2res)", "O2 -> O2(v3)", "O2 -> O2(v4)", "O2 -> O2(b1)",
            "O2 -> O2(4.5eV)", "O2 -> O2(9.97eV)"):
    i = proc_index("EXCITATION", tgt)
    e_tab(tgt.split("(")[1].rstrip(")"), "electron + O2 -> electron + O2",
          bz("EXCITATION", tgt), [p[3] for p in proc if p[0] == i][0], BZ + "Phelps O2, " + tgt)

# Heavy-particle reactions of Kimura, association rates at the operating p
PRs = "Plumb and Ryan (1986), Tables I-II, at the operating pressure"
heavy("k6", "CF3 + F -> CF4", pr(6), PRs + " (replaces k6)", rate=f"{pr(6):.3g} (fall-off)")
heavy("k7", "CF2 + F -> CF3", pr(8), PRs + " (replaces k7)", rate=f"{pr(8):.3g} (fall-off)")
heavy("k11", "CF + F -> CF2", pr(9), PRs + " (replaces k11)", rate=f"{pr(9):.3g} (fall-off)")
heavy("k17", "CF3 + Fn1 -> CF4 + electron", 5e-16, KN)
for x in POS:
    heavy("k18", f"{x} + Fn1 -> {NEUTRAL_OF[x]} + F", 4e-13, KN + " (X+ + F-)")
heavy("k24", "2F -> F2", 6.77e-40 * M, KN + " (three-body, times M)", order=2,
      rate=f"6.77e-40 M = {6.77e-40 * M:.3g}")
heavy("k27", "Fp1 + CF4 -> CF3p1 + F2", 1e-15, KN)
heavy("k28", "Fp1 + CF3 -> CF3p1 + F", 1e-15, KN)
heavy("k32", "On1 + Op1 -> 2O", 2.7e-13, KN)
for x in POS:
    if x != "Op1":
        heavy("k33", f"{x} + On1 -> {NEUTRAL_OF[x]} + O", 1.5e-13, KN + " (X+ + O-)")
heavy("k40", "On1 + O -> O2 + electron", 3e-16, KN)
heavy("k46", "On1 + O2a -> O2n1 + O", 1.0e-17, KN)
heavy("k47", "On1 + O2a -> O3 + electron", 2e-17, KN)
heavy("k48", "O + CF3 -> COF2 + F", 3.1e-17, KN)
heavy("k49", "O + CF2 -> COF + F", 1.4e-17, KN)
heavy("k50", "O + CF -> CO + F", 2.0e-17, KN)
heavy("k51", "O + COF -> CO2 + F", 9.3e-17, KN)
heavy("k52", "O + FO -> O2 + F", 5.0e-17, KN)
heavy("k53", "F + O -> FO", 1.0e-39 * M, KN + " (three-body, times M)", order=2,
      rate=f"1.0e-39 M = {1.0e-39 * M:.3g}")
heavy("k54", "CF2 + COF -> CF3 + CO", 3.0e-19, KN)
heavy("k55", "CF2 + COF -> COF2 + CF", 3.0e-19, KN)
heavy("k56", "CF3 + COF -> CF4 + CO", 1.0e-17, KN)
heavy("k57", "CF3 + COF -> COF2 + CF2", 1.0e-17, KN)
heavy("k58", "COF + COF -> COF2 + CO", 1.0e-17, KN)
heavy("k64", "Op1 + CF4 -> CF3p1 + FO", 1.4e-15, KN)
# Plumb and Ryan reactions that Kimura does not have
heavy("PR29", "CF2 + O -> CO + 2F", 4e-14 * 1e-6, "Plumb and Ryan (1986), Table I, reaction 29")
heavy("PR31", "COF + F -> COF2", pr(31), PRs + ", reaction 31", rate=f"{pr(31):.3g} (fall-off)")
heavy("PR38", "F + CO -> COF", pr(38), PRs + ", reaction 38", rate=f"{pr(38):.3g} (fall-off)")
# Oxygen chemistry of Toneli et al. (2015) and Atkinson et al. (2004)
TO = "Toneli et al. (2015), reaction "
k_o3 = 6.0e-34 * 1e-12 * (TG / 300.0) ** -2.6 * M
heavy("IUPAC", "O + O2 -> O3", k_o3, "Atkinson et al. (2004), O + O2 + M, M = O2, times M", order=2,
      rate=f"6.0e-46 (Tg/300)^-2.6 M = {k_o3:.3g}")
heavy("T64", "O2 + O2a -> 2O2", 2.2e-24, TO + "64", b=0.8)
heavy("T73", "O2a + O3 -> O + 2O2", 6.01e-17, TO + "73", Ta=2853)
heavy("T77", "O2a + O2a -> 2O2", 5.5e-29, TO + "77", b=0.5)
heavy("T82", "O + O3 -> 2O2", 1.81e-17, TO + "82", Ta=2300)
heavy("T84", "O + O3 -> O2 + O2a", 1e-17, TO + "84", Ta=2300)
heavy("T86", "O + O2a -> O + O2", 1.3e-22, TO + "86")
heavy("T100", "O2n1 + O2 -> electron + 2O2", 2.7e-16, TO + "100", b=0.5, Ta=5590)
heavy("T101", "O2n1 + O2a -> electron + 2O2", 2e-16, TO + "101")
heavy("T103", "O2n1 + O -> O2 + On1", 1.755e-16, TO + "103")
heavy("T104", "O2n1 + O -> electron + O3", 2.145e-16, TO + "104")
heavy("T105", "O2n1 + O2p1 -> 2O2", 2e-14, TO + "105", b=0.5)
heavy("T106", "O2n1 + O2p1 -> 2O + O2", 2e-14, TO + "106", b=0.5)
heavy("T107", "O2n1 + Op1 -> O + O2", 2e-14, TO + "107", b=0.5)
for x in POS:
    if x not in ("O2p1", "Op1"):
        heavy("T105+106", f"{x} + O2n1 -> {NEUTRAL_OF[x]} + O2", 4e-14,
              TO + "105 + 106, applied to X+ + O2-", b=0.5)
heavy("T109", "On1 + O3 -> electron + 2O2", 3.01e-16, TO + "109", b=0.5)
heavy("T110", "On1 + O3 -> O2 + O2n1", 0.3e-15, TO + "110")
heavy("T121", "O2 + Op1 -> O + O2p1", 2e-17, TO + "121", b=-0.5)
heavy("T122", "O3 + Op1 -> O2 + O2p1", 1e-16, TO + "122")
heavy("T123", "O2a + Op1 -> O + O2p1", 2e-17, TO + "123")
heavy("T126", "O2a -> O2", 2.256e-4, TO + "126 (radiative)")

# --------------------------------------------------------------------------
# Element and charge balance
for r in R:
    for key in ("C", "F", "O"):
        a = sum(COMP[s].get(key, 0) for s in r["lhs"])
        b = sum(COMP[s].get(key, 0) for s in r["rhs"])
        assert a == b, (r["label"], key)
    assert sum(CHARGE[s] for s in r["lhs"]) == sum(CHARGE[s] for s in r["rhs"]), r["label"]
print(f"{len(NAMES)} species ({N_CHARGED} charged), {len(R)} reactions, "
      f"M = {M:.5g} m^-3; elements and charge balance")


# --------------------------------------------------------------------------
# Output
def eqn(r):
    return " + ".join(r["lhs"]) + " = " + " + ".join(r["rhs"])


HDR = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\\\    /   O peration     | Version:     4.0                                |
|   \\\\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      %s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""
out = [HDR % "reactions",
       f"// CF4/O2 mechanism written by build_cf4o2_mechanism.py for p = {P:g} Pa\n"
       f"// and Tg = {TG:g} K; sources in MECHANISM.md\n\n",
       "species\n(\n" + "".join(f"    {s}\n" for s in NAMES) + ");\n\nreactions\n(\n"]
for r in R:
    out.append(f"    // {r['label']}: {r['src']}\n    {r['type']}\n    {eqn(r)}\n")
    if "table" in r:
        out.append(f'    ("$FOAM_CASE/constant/{r["table"]}")\n')
    else:
        out.append(f"    ({r['A']:.10e} {r['b']:.6g} {r['Ta']:.10g})\n")
    if r["type"] in ("electronImpactInelasticArrheniusReaction",
                     "electronImpactInelasticTabularReaction"):
        out.append(f"    ({r['dE']:.10e})\n")
    out.append("\n")
out.append(");\n")
(const / "reactions").write_text("".join(out))

for name, t in tables.items():
    (const / name).write_text("(\n" + "".join(f"({a:.10g} {b:.10g})\n" for a, b in t) + ")\n")
for name in ("mu_electron", "D_electron"):
    t = read_table(bdir / name)
    (const / name).write_text("(\n" + "".join(f"{a:.10g} {b:.10g}\n" for a, b in t) + ")\n")

th = ["(\n"]
for s in NAMES:
    th.append(f"{s} {s} 1 {molar_mass(s):.8g}\n    20.79 1.188e+06\n    2.09e-05 0.7\n")
th.append(")\n")
(const / "thermoProperties").write_text("".join(th))

md = ["# CF4/O2 mechanism\n\n",
      f"Written by `build_cf4o2_mechanism.py` for p = {P:g} Pa and Tg = {TG:g} K "
      f"(M = {M:.5g} m^-3). {len(NAMES)} species, {N_CHARGED} of them charged, "
      f"and {len(R)} reactions. Rates in m3/s (1/s for first order), Te in eV, "
      "dE the electron energy loss per event in eV. The choices are described in the "
      "header of the script.\n\n",
      "Species: " + ", ".join(NAMES) + ".\n\n",
      "| # | Label | Reaction | Rate | dE (eV) | Source |\n|---|---|---|---|---|---|\n"]
for i, r in enumerate(R, 1):
    dE = "" if r["dEeV"] is None else ("elastic" if "Elastic" in r["type"] else f"{r['dEeV']:g}")
    md.append(f"| {i} | {r['label']} | {' + '.join(r['lhs'])} → {' + '.join(r['rhs'])} "
              f"| {r['rate']} | {dE} | {r['src']} |\n")
(case / "MECHANISM.md").write_text("".join(md))
print("written", const / "reactions", "and", case / "MECHANISM.md")
