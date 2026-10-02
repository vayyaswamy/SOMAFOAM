# SOMAFOAM

This software consists of FOAM based FV method and other utilities designed and implemented for modular, multiphysics plasma fluid simulation.

The foam-extend snapshot used in this project corresponds to following,
```
commit efcc2b1b7df8543c7873f89a6f50c30e047b1b11
Date:   Thu Apr 12 14:05:28 2018 +0100
```
There were certain changes made to the foam-extend base depending on our needs.
The code compiles on newer gcc versions (including 10.2.1)

## Building

The repository contains sources only; compiled libraries and applications are
machine-specific (compiler, glibc) and are built locally.

Requirements: gcc/g++, make, flex and OpenMPI (with development headers).
Tested on Debian 11 (gcc 10.2.1, OpenMPI 4.1.0).

```
git clone https://github.com/vayyaswamy/SOMAFOAM
cd SOMAFOAM
./install.sh            # builds wmake, all libraries and applications
```

The libraries go to `lib/linux64GccDPOpt` and the applications to
`bin/linux64GccDPOpt`. Before running, load the environment in each new shell:

```
source <path-to>/SOMAFOAM/etc/bashrc
```

Run the cases in `examples/plasma` with `somaFoam` (their `controlDict` still
names the older `plasmaSimFoam`, which no longer runs them).

Contributors:
1) Venkattraman Ayyaswamy (https://me.ucmerced.edu/content/venkattraman-venkatt-ayyaswamy)
2) Abhishek Kumar Verma 
3) Saurav Gautam (https://me.ucmerced.edu/content/saurav-gautam)
4) Jose Alfredo Millan Higuera (https://www.linkedin.com/in/jmillanhiguera/)

For research conducted using SOMAFOAM please cite this paper: https://doi.org/10.1016/j.cpc.2021.107855

Paper that compares SOMAFOAM results with experiment: 
https://doi.org/10.1063/5.0041386
