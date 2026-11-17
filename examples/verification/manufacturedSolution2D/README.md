# Manufactured solution in 2D (paper case V1)

A steady manufactured-solution test of the discretisation that `somaFoam` uses
for the potential, the drift-diffusion of the ions and the electron energy.
The test application `mmsPlasma2D` (`applications/verification/mmsPlasma2D`)
solves, on the unit square,

    -laplacian(V)                     = S_V
    div(n u_i) - laplacian(D_i, n)    = S_n,  u_i = -mu_i grad(V)
    div(w u_e) - laplacian(D_e, w)    = S_w,  u_e =  mu_e grad(V)

with D_i = 1, the variable coefficient D_e = 1 + x y/2, mu_i = 20 and mu_e = 40
(`system/mmsDict`). It builds the drift fluxes and the matrices with the same
functions as `multiSpeciesPlasmaModel`, either the Scharfetter-Gummel flux or
the `fvSchemes` discretisation with `Gauss limitedLinear 1` convection and the
corrected Laplacian of the plasma examples. The exact fields are

    V = sin(pi x) cos(pi y)/2 + x y
    n = exp(x/2) (1 + 0.4 sin(pi x) sin(pi y))
    w = 1.5 + cos(pi x) cos(pi y)/2 + y/5

and `generate.py` (sympy) writes them and the sources to
`applications/verification/mmsPlasma2D/mmsSolution.H`. The boundary values are
the exact solution.

`./Allrun` runs 16, 32, 64 and 128 cells per side on a uniform mesh and on a
mesh distorted by `distortMesh.py` (interior points moved by
0.06 sin(2 pi x) sin(2 pi y) in x and y, maximum non-orthogonality 31 degrees,
skewness 0.38). `verify.py` prints the L2 errors and the observed orders,
writes `runs/mmsOrder.pdf`, and prints PASS when the orders on the two finest
meshes reach the expected orders within 0.2. The whole run took 11 s on
one core.

Results of 2026-11-17 on Pinnacles, L2 error and order between 64 and 128
cells per side:

| Mesh | Flux | V | n | w |
|---|---|---|---|---|
| uniform | Scharfetter-Gummel | 1.69e-5, 2.00 | 5.02e-4, 1.94 | 8.97e-4, 1.88 |
| uniform | fvSchemes | 1.69e-5, 2.00 | 3.88e-4, 1.95 | 1.07e-3, 1.86 |
| distorted | Scharfetter-Gummel | 7.50e-4, 0.99 | 3.80e-2, 0.70 | 2.41e-2, 1.08 |
| distorted | fvSchemes | 7.50e-4, 0.99 | 2.82e-2, 0.75 | 3.32e-2, 1.04 |

The cell Peclet number of the drift fell from 3 to 0.4 for the ions and from
4.3 to 0.5 for the electrons over the four meshes, so the coarse meshes are
drift dominated and the fine ones diffusion dominated. On the distorted mesh
the solution was first order with both flux schemes. The order of the
potential stayed at 1.0 with the least-squares gradient and with 20 instead of
4 non-orthogonal correctors, which places the first-order error in the face
gradient of the non-orthogonal correction and in the linear interpolation
without skewness correction rather than in the iteration. The order of the
ion density was still rising (0.4, 0.6, 0.75) and is reported without a check.
