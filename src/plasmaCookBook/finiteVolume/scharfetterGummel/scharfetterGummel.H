/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Function
    Foam::scharfetterGummel

Description
    Implicit matrix of the divergence of a drift-diffusion flux,

        div(phi psi) - laplacian(gamma, psi),

    discretised with the Scharfetter-Gummel (exponential) scheme: drift and
    diffusion are treated as one flux, which between two neighbouring cell
    centres P and N is the exact solution of the steady one-dimensional
    problem,

        flux = gamma |Sf|/d [B(-Pe) psi_P - B(Pe) psi_N],

        Pe = phi d/(gamma |Sf|),    B(x) = x/(exp(x) - 1).

    It reduces to central differencing for small Pe and to upwind
    differencing for large |Pe|, and gives non-negative solutions at any
    cell Peclet number without a limiter.

    phi is the volumetric drift flux through the faces and gamma the
    diffusion coefficient, interpolated linearly to the faces. The two-point
    formula is used on internal faces and on processor and cyclic patches;
    other boundary faces are treated as by the Gauss convection and
    Laplacian schemes, so that the boundary conditions act as before. There
    is no non-orthogonal correction.

SourceFiles
    scharfetterGummel.C

\*---------------------------------------------------------------------------*/

#ifndef scharfetterGummel_H
#define scharfetterGummel_H

#include "fvMatrices.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<fvScalarMatrix> scharfetterGummel
(
    const surfaceScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& psi
);

//- The flux through the faces that the matrix above stands for, evaluated
//  for a given psi
tmp<surfaceScalarField> scharfetterGummelFlux
(
    const surfaceScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& psi
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
