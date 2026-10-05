/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "scharfetterGummel.H"
#include "linear.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Coefficients of the flux out of the owner cell,
//  flux = aOwn psiOwn - aNei psiNei, with aOwn = g B(-Pe), aNei = g B(Pe),
//  Pe = phi/g and g = gamma |Sf|/d
static inline void scharfetterGummelCoeffs
(
    const scalar gIn,
    const scalar phi,
    scalar& aOwn,
    scalar& aNei
)
{
    const scalar g = max(gIn, scalar(0));

    if (mag(phi) <= 1e-4*g)
    {
        // Diffusion dominated: series of the Bernoulli function
        aOwn = g + 0.5*phi;
        aNei = g - 0.5*phi;
    }
    else if (mag(phi) >= 100*g)
    {
        // Drift dominated (or no diffusion): upwind
        aOwn = max(phi, scalar(0));
        aNei = max(-phi, scalar(0));
    }
    else
    {
        aNei = phi/(Foam::exp(phi/g) - 1.0);
        aOwn = aNei + phi;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::scharfetterGummel
(
    const surfaceScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& psi
)
{
    const fvMesh& mesh = psi.mesh();

    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix(psi, phi.dimensions()*psi.dimensions())
    );
    fvScalarMatrix& fvm = tfvm();

    // gamma |Sf|/d on the faces
    const surfaceScalarField g
    (
        linear<scalar>(mesh).interpolate(gamma)*mesh.magSf()*mesh.deltaCoeffs()
    );

    scalarField& lower = fvm.lower();
    scalarField& upper = fvm.upper();

    const scalarField& gI = g.internalField();
    const scalarField& phiI = phi.internalField();

    forAll(gI, faceI)
    {
        scalar aOwn, aNei;
        scharfetterGummelCoeffs(gI[faceI], phiI[faceI], aOwn, aNei);

        lower[faceI] = -aOwn;
        upper[faceI] = -aNei;
    }

    fvm.negSumDiag();

    forAll(psi.boundaryField(), patchI)
    {
        const fvPatchScalarField& ppsi = psi.boundaryField()[patchI];
        const scalarField& pphi = phi.boundaryField()[patchI];

        scalarField& internalCoeffs = fvm.internalCoeffs()[patchI];
        scalarField& boundaryCoeffs = fvm.boundaryCoeffs()[patchI];

        if (ppsi.coupled() && ppsi.patch().type() != "regionCouple")
        {
            const scalarField& pg = g.boundaryField()[patchI];

            forAll(pg, faceI)
            {
                scharfetterGummelCoeffs
                (
                    pg[faceI],
                    pphi[faceI],
                    internalCoeffs[faceI],
                    boundaryCoeffs[faceI]
                );
            }
        }
        else
        {
            const fvsPatchScalarField& pw = mesh.weights().boundaryField()[patchI];

            const scalarField pGammaMagSf
            (
                gamma.boundaryField()[patchI]*ppsi.patch().magSf()
            );

            internalCoeffs =
                pphi*ppsi.valueInternalCoeffs(pw)
              - pGammaMagSf*ppsi.gradientInternalCoeffs();

            boundaryCoeffs =
              - pphi*ppsi.valueBoundaryCoeffs(pw)
              + pGammaMagSf*ppsi.gradientBoundaryCoeffs();
        }
    }

    return tfvm;
}


Foam::tmp<Foam::surfaceScalarField> Foam::scharfetterGummelFlux
(
    const surfaceScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& psi
)
{
    const fvMesh& mesh = psi.mesh();

    tmp<surfaceScalarField> tflux
    (
        new surfaceScalarField
        (
            IOobject
            (
                "scharfetterGummelFlux(" + psi.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", phi.dimensions()*psi.dimensions(), 0.0)
        )
    );
    surfaceScalarField& flux = tflux();

    const surfaceScalarField g
    (
        linear<scalar>(mesh).interpolate(gamma)*mesh.magSf()*mesh.deltaCoeffs()
    );

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const scalarField& gI = g.internalField();
    const scalarField& phiI = phi.internalField();
    const scalarField& psiI = psi.internalField();
    scalarField& fluxI = flux.internalField();

    forAll(fluxI, faceI)
    {
        scalar aOwn, aNei;
        scharfetterGummelCoeffs(gI[faceI], phiI[faceI], aOwn, aNei);

        fluxI[faceI] = aOwn*psiI[owner[faceI]] - aNei*psiI[neighbour[faceI]];
    }

    forAll(psi.boundaryField(), patchI)
    {
        const fvPatchScalarField& ppsi = psi.boundaryField()[patchI];
        const scalarField& pphi = phi.boundaryField()[patchI];
        scalarField& pflux = flux.boundaryField()[patchI];

        if (ppsi.coupled() && ppsi.patch().type() != "regionCouple")
        {
            const scalarField& pg = g.boundaryField()[patchI];
            const scalarField psiOwn(ppsi.patchInternalField());
            const scalarField psiNei(ppsi.patchNeighbourField());

            forAll(pflux, faceI)
            {
                scalar aOwn, aNei;
                scharfetterGummelCoeffs(pg[faceI], pphi[faceI], aOwn, aNei);

                pflux[faceI] = aOwn*psiOwn[faceI] - aNei*psiNei[faceI];
            }
        }
        else if (pflux.size())
        {
            pflux =
                pphi*ppsi
              - gamma.boundaryField()[patchI]*ppsi.patch().magSf()
               *ppsi.snGrad();
        }
    }

    return tflux;
}


// ************************************************************************* //
