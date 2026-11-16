# CF4/O2 mechanism

Written by `build_cf4o2_mechanism.py` for p = 133.322 Pa and Tg = 300 K (M = 3.2188e+22 m^-3). 29 species, 13 of them charged, and 125 reactions. Rates in m3/s (1/s for first order), Te in eV, dE the electron energy loss per event in eV. The choices are described in the header of the script.

Species: electron, CF3p1, CF2p1, CFp1, Fp1, O2p1, Op1, COp1, CO2p1, Cp1, Fn1, On1, O2n1, CF3, CF2, CF, F, F2, O2, O2a, O, O3, CO, CO2, FO, COF, COF2, C, CF4.

| # | Label | Reaction | Rate | dE (eV) | Source |
|---|---|---|---|---|---|
| 1 | k1 | electron + CF4 → electron + CF3 + F | cross section over the EEDF | 12.85 | LXCat Hayashi CF4 neutral dissociation (12.85 eV) times the CF3 fraction of Christophorou et al. (1996), Table 21 (replaces Kimura k1) |
| 2 | k2 | electron + CF4 → electron + CF2 + F + F | cross section over the EEDF | 12.85 | LXCat Hayashi CF4 neutral dissociation (12.85 eV) times the CF2 fraction of Christophorou et al. (1996), Table 21 (replaces Kimura k2) |
| 3 | k3 | electron + CF4 → electron + CF + F + F + F | cross section over the EEDF | 12.85 | LXCat Hayashi CF4 neutral dissociation (12.85 eV) times the CF1 fraction of Christophorou et al. (1996), Table 21 (replaces Kimura k3) |
| 4 | k8 | electron + CF4 → electron + electron + CF3p1 + F | cross section over the EEDF | 16.18 | LXCat Hayashi CF4 total ionization times the CF3+ fraction of Christophorou et al. (1996), Table 15 (replaces Kimura k8) |
| 5 | k9 | electron + CF4 → electron + electron + CF2p1 + F + F | cross section over the EEDF | 16.18 | LXCat Hayashi CF4 total ionization times the CF2+ fraction of Christophorou et al. (1996), Table 15 (replaces Kimura k9) |
| 6 | k10 | electron + CF4 → electron + electron + CFp1 + F + F + F | cross section over the EEDF | 16.18 | LXCat Hayashi CF4 total ionization times the CF+ fraction of Christophorou et al. (1996), Table 15 (replaces Kimura k10) |
| 7 | CXi | electron + CF4 → electron + electron + Cp1 + F + F + F + F | cross section over the EEDF | 16.18 | LXCat Hayashi CF4 total ionization times the C+ fraction of Christophorou et al. (1996), Table 15 |
| 8 | FXi | electron + CF4 → electron + electron + Fp1 + CF3 | cross section over the EEDF | 16.18 | LXCat Hayashi CF4 total ionization times the F+ fraction of Christophorou et al. (1996), Table 15 |
| 9 | k15 | electron + CF4 → Fn1 + CF3 | cross section over the EEDF | 0 | LXCat Hayashi CF4 total attachment (F- and CF3-) as F- + CF3 (replaces Kimura k15) |
| 10 | k30 | electron + O2 → electron + electron + O2p1 | cross section over the EEDF | 12.06 | LXCat Phelps O2 O2 -> O2^+ (replaces Kimura k30) |
| 11 | k31 | electron + O2 → O + On1 | cross section over the EEDF | 0 | LXCat Phelps O2 O2 -> O^-+O (replaces Kimura k31) |
| 12 | k35 | electron + O2 → electron + O + O | cross section over the EEDF | 6 | LXCat Phelps O2 O2 -> O2(6.0eV) (replaces Kimura k35) |
| 13 | k36 | electron + O2 → electron + O + O | cross section over the EEDF | 8.4 | LXCat Phelps O2 O2 -> O2(8.4eV) (replaces Kimura k36) |
| 14 | k41 | electron + O2 → electron + O2a | cross section over the EEDF | 0.977 | LXCat Phelps O2 O2 -> O2(a1) (replaces Kimura k41) |
| 15 | k44 | electron + O2a → electron + O2 | cross section over the EEDF | -0.977 | LXCat Phelps O2 O2(a1) by detailed balance (replaces Kimura k44) |
| 16 | k4 | electron + CF3 → electron + CF2 + F | 6.48e-14 Te^-0.96 exp(-11.25/Te) | 11.25 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 17 | k5 | electron + CF2 → electron + CF + F | 1.16e-14 Te^-0.38 exp(-14.35/Te) | 14.35 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 18 | k12 | electron + CF3 → electron + electron + CF3p1 | 8e-15 exp(-12.2/Te) | 12.2 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 19 | k13 | electron + CF2 → electron + electron + CF2p1 | 2.5e-14 exp(-12.2/Te) | 12.2 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 20 | k14 | electron + CF → electron + electron + CFp1 | 2.5e-14 exp(-15.3/Te) | 15.3 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 21 | k16a | electron + CF2p1 → CF2 | 6e-14 | 0 | Kimura and Noto (2006), Table I |
| 22 | k16b | electron + CF3p1 → CF3 | 6e-14 | 0 | Kimura and Noto (2006), Table I |
| 23 | k19 | electron + CF3 → electron + electron + CF2p1 + F | 2.2e-14 exp(-18.6/Te) | 18.6 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 24 | k20 | electron + CF2 → electron + electron + CFp1 + F | 2.33e-14 exp(-17/Te) | 17 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 25 | k21 | electron + CF3 → CF2 + Fn1 | 1.6e-16 Te^-0.4329 exp(-0.18/Te) | 0 | Kimura and Noto (2006), Table I |
| 26 | k22 | electron + CF2 → CF + Fn1 | 1.6e-16 Te^-0.4329 exp(-1.93/Te) | 0 | Kimura and Noto (2006), Table I |
| 27 | k23 | electron + F → electron + electron + Fp1 | 1.3e-14 exp(-16.5/Te) | 16.5 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 28 | k25 | electron + F2 → Fn1 + F | 4.5e-15 Te^-1.35 exp(-0.15/Te) | 0 | Kimura and Noto (2006), Table I |
| 29 | k26 | electron + F2 → electron + F + F | 1.18e-14 exp(-5.77/Te) | 5.77 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 30 | k29 | electron + O → electron + electron + Op1 | 9e-15 Te^0.7 exp(-13.6/Te) | 13.6 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 31 | k34 | electron + On1 → electron + electron + O | 2e-13 exp(-5.5/Te) | 5.5 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 32 | k39 | electron + O2p1 → O + O | 5.2e-15 Te^-1 exp(-0/Te) | 0 | Kimura and Noto (2006), Table I |
| 33 | k42 | electron + O2a → electron + electron + O2p1 | 9e-16 Te^2 exp(-11.6/Te) | 11.6 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 34 | k43 | electron + O2a → On1 + O | 2.28e-16 exp(-2.29/Te) | 0 | Kimura and Noto (2006), Table I |
| 35 | k45 | electron + O2a → electron + O + O | 4.2e-15 Te^2 exp(-4.6/Te) | 4.6 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 36 | k59 | electron + COF2 → electron + COF + F | 6.48e-14 Te^-0.96 exp(-11.25/Te) | 11.25 | Kimura and Noto (2006), Table I (k59 = k4); dE activation energy of the fit |
| 37 | k60 | electron + CO2 → electron + CO + O | 8.5e-16 | 0 | Kimura and Noto (2006), Table I (constant) |
| 38 | k61 | electron + CO → electron + C + O | 1.87e-14 exp(-13.89/Te) | 13.89 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 39 | k62 | electron + CO2 → electron + electron + CO2p1 | 5.03e-14 exp(-16.26/Te) | 16.26 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 40 | k63 | electron + CO → electron + electron + COp1 | 4.93e-14 exp(-16.76/Te) | 16.76 | Kimura and Noto (2006), Table I; dE activation energy of the fit |
| 41 | elCF4 | electron + CF4 → electron + CF4 | Boltzmann table | elastic | Boltzmann equation, LXCat Hayashi CF4 |
| 42 | elO2 | electron + O2 → electron + O2 | Boltzmann table | elastic | Boltzmann equation, LXCat Phelps O2, effective minus inelastic |
| 43 | V24 | electron + CF4 → electron + CF4 | Boltzmann table | 0.052 | Boltzmann equation, LXCat Hayashi CF4, CF4 -> CF4(V24) |
| 44 | V13 | electron + CF4 → electron + CF4 | Boltzmann table | 0.109 | Boltzmann equation, LXCat Hayashi CF4, CF4 -> CF4(V13) |
| 45 | 7.86eV | electron + CF4 → electron + CF4 | Boltzmann table | 7.86 | Boltzmann equation, LXCat Hayashi CF4, CF4 -> CF4*(7.86eV) |
| 46 | rot | electron + O2 → electron + O2 | Boltzmann table | 0.02 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(rot) |
| 47 | v1 | electron + O2 → electron + O2 | Boltzmann table | 0.19 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(v1) |
| 48 | v1res | electron + O2 → electron + O2 | Boltzmann table | 0.19 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(v1res) |
| 49 | v2 | electron + O2 → electron + O2 | Boltzmann table | 0.38 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(v2) |
| 50 | v2res | electron + O2 → electron + O2 | Boltzmann table | 0.38 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(v2res) |
| 51 | v3 | electron + O2 → electron + O2 | Boltzmann table | 0.57 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(v3) |
| 52 | v4 | electron + O2 → electron + O2 | Boltzmann table | 0.75 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(v4) |
| 53 | b1 | electron + O2 → electron + O2 | Boltzmann table | 1.627 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(b1) |
| 54 | 4.5eV | electron + O2 → electron + O2 | Boltzmann table | 4.5 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(4.5eV) |
| 55 | 9.97eV | electron + O2 → electron + O2 | Boltzmann table | 9.97 | Boltzmann equation, LXCat Phelps O2, O2 -> O2(9.97eV) |
| 56 | k6 | CF3 + F → CF4 | 1.5e-17 (fall-off) |  | Plumb and Ryan (1986), Tables I-II, at the operating pressure (replaces k6) |
| 57 | k7 | CF2 + F → CF3 | 7.83e-19 (fall-off) |  | Plumb and Ryan (1986), Tables I-II, at the operating pressure (replaces k7) |
| 58 | k11 | CF + F → CF2 | 9.95e-21 (fall-off) |  | Plumb and Ryan (1986), Tables I-II, at the operating pressure (replaces k11) |
| 59 | k17 | CF3 + Fn1 → CF4 + electron | 5e-16 |  | Kimura and Noto (2006), Table I |
| 60 | k18 | CF3p1 + Fn1 → CF3 + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 61 | k18 | CF2p1 + Fn1 → CF2 + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 62 | k18 | CFp1 + Fn1 → CF + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 63 | k18 | Fp1 + Fn1 → F + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 64 | k18 | O2p1 + Fn1 → O2 + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 65 | k18 | Op1 + Fn1 → O + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 66 | k18 | COp1 + Fn1 → CO + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 67 | k18 | CO2p1 + Fn1 → CO2 + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 68 | k18 | Cp1 + Fn1 → C + F | 4e-13 |  | Kimura and Noto (2006), Table I (X+ + F-) |
| 69 | k24 | F + F → F2 | 6.77e-40 M = 2.18e-17 |  | Kimura and Noto (2006), Table I (three-body, times M) |
| 70 | k27 | Fp1 + CF4 → CF3p1 + F2 | 1e-15 |  | Kimura and Noto (2006), Table I |
| 71 | k28 | Fp1 + CF3 → CF3p1 + F | 1e-15 |  | Kimura and Noto (2006), Table I |
| 72 | k32 | On1 + Op1 → O + O | 2.7e-13 |  | Kimura and Noto (2006), Table I |
| 73 | k33 | CF3p1 + On1 → CF3 + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 74 | k33 | CF2p1 + On1 → CF2 + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 75 | k33 | CFp1 + On1 → CF + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 76 | k33 | Fp1 + On1 → F + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 77 | k33 | O2p1 + On1 → O2 + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 78 | k33 | COp1 + On1 → CO + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 79 | k33 | CO2p1 + On1 → CO2 + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 80 | k33 | Cp1 + On1 → C + O | 1.5e-13 |  | Kimura and Noto (2006), Table I (X+ + O-) |
| 81 | k40 | On1 + O → O2 + electron | 3e-16 |  | Kimura and Noto (2006), Table I |
| 82 | k46 | On1 + O2a → O2n1 + O | 1e-17 |  | Kimura and Noto (2006), Table I |
| 83 | k47 | On1 + O2a → O3 + electron | 2e-17 |  | Kimura and Noto (2006), Table I |
| 84 | k48 | O + CF3 → COF2 + F | 3.1e-17 |  | Kimura and Noto (2006), Table I |
| 85 | k49 | O + CF2 → COF + F | 1.4e-17 |  | Kimura and Noto (2006), Table I |
| 86 | k50 | O + CF → CO + F | 2e-17 |  | Kimura and Noto (2006), Table I |
| 87 | k51 | O + COF → CO2 + F | 9.3e-17 |  | Kimura and Noto (2006), Table I |
| 88 | k52 | O + FO → O2 + F | 5e-17 |  | Kimura and Noto (2006), Table I |
| 89 | k53 | F + O → FO | 1.0e-39 M = 3.22e-17 |  | Kimura and Noto (2006), Table I (three-body, times M) |
| 90 | k54 | CF2 + COF → CF3 + CO | 3e-19 |  | Kimura and Noto (2006), Table I |
| 91 | k55 | CF2 + COF → COF2 + CF | 3e-19 |  | Kimura and Noto (2006), Table I |
| 92 | k56 | CF3 + COF → CF4 + CO | 1e-17 |  | Kimura and Noto (2006), Table I |
| 93 | k57 | CF3 + COF → COF2 + CF2 | 1e-17 |  | Kimura and Noto (2006), Table I |
| 94 | k58 | COF + COF → COF2 + CO | 1e-17 |  | Kimura and Noto (2006), Table I |
| 95 | k64 | Op1 + CF4 → CF3p1 + FO | 1.4e-15 |  | Kimura and Noto (2006), Table I |
| 96 | PR29 | CF2 + O → CO + F + F | 4e-20 |  | Plumb and Ryan (1986), Table I, reaction 29 |
| 97 | PR31 | COF + F → COF2 | 1.45e-18 (fall-off) |  | Plumb and Ryan (1986), Tables I-II, at the operating pressure, reaction 31 |
| 98 | PR38 | F + CO → COF | 2.57e-21 (fall-off) |  | Plumb and Ryan (1986), Tables I-II, at the operating pressure, reaction 38 |
| 99 | IUPAC | O + O2 → O3 | 6.0e-46 (Tg/300)^-2.6 M = 1.93e-23 |  | Atkinson et al. (2004), O + O2 + M, M = O2, times M |
| 100 | T64 | O2 + O2a → O2 + O2 | 2.2e-24 (Tg/300)^0.8 |  | Toneli et al. (2015), reaction 64 |
| 101 | T73 | O2a + O3 → O + O2 + O2 | 6.01e-17 exp(-2853/Tg) |  | Toneli et al. (2015), reaction 73 |
| 102 | T77 | O2a + O2a → O2 + O2 | 5.5e-29 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 77 |
| 103 | T82 | O + O3 → O2 + O2 | 1.81e-17 exp(-2300/Tg) |  | Toneli et al. (2015), reaction 82 |
| 104 | T84 | O + O3 → O2 + O2a | 1e-17 exp(-2300/Tg) |  | Toneli et al. (2015), reaction 84 |
| 105 | T86 | O + O2a → O + O2 | 1.3e-22 |  | Toneli et al. (2015), reaction 86 |
| 106 | T100 | O2n1 + O2 → electron + O2 + O2 | 2.7e-16 (Tg/300)^0.5 exp(-5590/Tg) |  | Toneli et al. (2015), reaction 100 |
| 107 | T101 | O2n1 + O2a → electron + O2 + O2 | 2e-16 |  | Toneli et al. (2015), reaction 101 |
| 108 | T103 | O2n1 + O → O2 + On1 | 1.75e-16 |  | Toneli et al. (2015), reaction 103 |
| 109 | T104 | O2n1 + O → electron + O3 | 2.15e-16 |  | Toneli et al. (2015), reaction 104 |
| 110 | T105 | O2n1 + O2p1 → O2 + O2 | 2e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 |
| 111 | T106 | O2n1 + O2p1 → O + O + O2 | 2e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 106 |
| 112 | T107 | O2n1 + Op1 → O + O2 | 2e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 107 |
| 113 | T105+106 | CF3p1 + O2n1 → CF3 + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 114 | T105+106 | CF2p1 + O2n1 → CF2 + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 115 | T105+106 | CFp1 + O2n1 → CF + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 116 | T105+106 | Fp1 + O2n1 → F + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 117 | T105+106 | COp1 + O2n1 → CO + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 118 | T105+106 | CO2p1 + O2n1 → CO2 + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 119 | T105+106 | Cp1 + O2n1 → C + O2 | 4e-14 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 105 + 106, applied to X+ + O2- |
| 120 | T109 | On1 + O3 → electron + O2 + O2 | 3.01e-16 (Tg/300)^0.5 |  | Toneli et al. (2015), reaction 109 |
| 121 | T110 | On1 + O3 → O2 + O2n1 | 3e-16 |  | Toneli et al. (2015), reaction 110 |
| 122 | T121 | O2 + Op1 → O + O2p1 | 2e-17 (Tg/300)^-0.5 |  | Toneli et al. (2015), reaction 121 |
| 123 | T122 | O3 + Op1 → O2 + O2p1 | 1e-16 |  | Toneli et al. (2015), reaction 122 |
| 124 | T123 | O2a + Op1 → O + O2p1 | 2e-17 |  | Toneli et al. (2015), reaction 123 |
| 125 | T126 | O2a → O2 | 0.000226 |  | Toneli et al. (2015), reaction 126 (radiative) |
