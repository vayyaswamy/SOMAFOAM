/*---------------------------------------------------------------------------*\
License
    This file is part of SOMAFOAM, distributed under the GNU General Public
    License, version 3 or later.

Application
    mmsPlasma2D

Description
    Steady manufactured-solution test of the discretisation used by somaFoam
    for the potential, the drift-diffusion of a positive ion and the
    electron energy density (examples/verification/manufacturedSolution2D).

    The equations, on the unit square, are

        -laplacian(V)                     = S_V
        div(n u_i) - laplacian(D_i, n)    = S_n,  u_i = -mu_i grad(V)
        div(w u_e) - laplacian(D_e, w)    = S_w,  u_e =  mu_e grad(V)

    with D_i = 1 and D_e = 1 + x y/2. The exact fields and the sources are in
    mmsSolution.H, written by generate.py of the example. The drift fluxes
    through the faces and the drift-diffusion matrices are built as in
    multiSpeciesPlasmaModel (driftFaceFlux and driftDiffusionTerms):

        driftFlux       cell: (interpolate(u) & Sf), with u from the cell
                        gradient of V; face: from snGrad(V)
        fluxScheme      scharfetterGummel, or fvSchemes for
                        div(phi,n) - laplacian(D,n) with the schemes of
                        system/fvSchemes

    Boundary values are the exact solution. The volume-weighted L2 and the
    maximum errors are printed and written to mmsErrors.dat.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "scharfetterGummel.H"
#include "OFstream.H"
#include "mmsSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

typedef double (*mmsFunction)(const double, const double);

void setCells(volScalarField& f, mmsFunction func)
{
    const vectorField& C = f.mesh().C().internalField();

    forAll(C, cellI)
    {
        f.internalField()[cellI] = func(C[cellI].x(), C[cellI].y());
    }
}

void setBoundary(volScalarField& f, mmsFunction func)
{
    forAll(f.boundaryField(), patchI)
    {
        fvPatchScalarField& pf = f.boundaryField()[patchI];

        if (pf.size() && !pf.coupled())
        {
            const vectorField& Cf = pf.patch().Cf();

            scalarField values(Cf.size());

            forAll(Cf, faceI)
            {
                values[faceI] = func(Cf[faceI].x(), Cf[faceI].y());
            }

            pf == values;
        }
    }
}

void setField(volScalarField& f, mmsFunction func)
{
    setCells(f, func);
    setBoundary(f, func);
}

void errors
(
    const volScalarField& f,
    mmsFunction exact,
    scalar& L2,
    scalar& Linf
)
{
    const vectorField& C = f.mesh().C().internalField();
    const scalarField& V = f.mesh().V();

    scalar sumE2 = 0;
    scalar sumV = 0;
    Linf = 0;

    forAll(C, cellI)
    {
        const scalar e = f[cellI] - exact(C[cellI].x(), C[cellI].y());

        sumE2 += V[cellI]*e*e;
        sumV += V[cellI];
        Linf = max(Linf, mag(e));
    }

    reduce(sumE2, sumOp<scalar>());
    reduce(sumV, sumOp<scalar>());
    reduce(Linf, maxOp<scalar>());

    L2 = Foam::sqrt(sumE2/sumV);
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary mmsDict
    (
        IOobject
        (
            "mmsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar muI(readScalar(mmsDict.lookup("muI")));
    const scalar muE(readScalar(mmsDict.lookup("muE")));
    const word fluxScheme(mmsDict.lookup("fluxScheme"));
    const word driftFlux(mmsDict.lookupOrDefault<word>("driftFlux", "cell"));
    const label nNonOrth
    (
        mmsDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0)
    );

    if (fluxScheme != "scharfetterGummel" && fluxScheme != "fvSchemes")
    {
        FatalErrorIn(args.executable())
            << "Unknown fluxScheme " << fluxScheme
            << "; valid entries are scharfetterGummel and fvSchemes"
            << exit(FatalError);
    }

    if (driftFlux != "cell" && driftFlux != "face")
    {
        FatalErrorIn(args.executable())
            << "Unknown driftFlux " << driftFlux
            << "; valid entries are cell and face" << exit(FatalError);
    }

    Info<< "fluxScheme " << fluxScheme << ", driftFlux " << driftFlux
        << ", mu_i " << muI << ", mu_e " << muE
        << ", nNonOrthogonalCorrectors " << nNonOrth << nl << endl;

    volScalarField V
    (
        IOobject("V", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField n
    (
        IOobject("n", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField w
    (
        IOobject("w", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );

    setBoundary(V, mms::Vexact);
    setBoundary(n, mms::nExact);
    setBoundary(w, mms::wExact);

    volScalarField Di
    (
        IOobject("Di", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Di", dimless, 1.0)
    );
    volScalarField De
    (
        IOobject("De", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("De", dimless, 0)
    );
    setField(De, mms::DeExact);

    // Sources at the cell centres, with the dimensions of the matrices
    volScalarField SV
    (
        IOobject("SV", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("SV", dimless/dimArea, 0)
    );
    setCells(SV, mms::SV);

    volScalarField Sn
    (
        IOobject("Sn", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Sn", dimless/dimArea, 0)
    );
    volScalarField Sw
    (
        IOobject("Sw", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Sw", dimless/dimArea, 0)
    );
    {
        const vectorField& C = mesh.C().internalField();

        forAll(C, cellI)
        {
            const scalar x = C[cellI].x();
            const scalar y = C[cellI].y();

            Sn[cellI] = mms::SnDiff(x, y) + muI*mms::SnDrift(x, y);
            Sw[cellI] = mms::SwDiff(x, y) + muE*mms::SwDrift(x, y);
        }
    }

    // Potential
    for (label corr = 0; corr <= nNonOrth; corr++)
    {
        fvScalarMatrix VEqn(-fvm::laplacian(V) == SV);
        VEqn.solve();
    }

    // Drift fluxes through the faces, as in driftFaceFlux
    surfaceScalarField phiI
    (
        IOobject("phiI", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("phiI", dimArea/dimLength, 0)
    );
    surfaceScalarField phiE(phiI);

    if (driftFlux == "cell")
    {
        const volVectorField gradV(fvc::grad(V));

        phiI = fvc::interpolate(-muI*gradV) & mesh.Sf();
        phiE = fvc::interpolate(muE*gradV) & mesh.Sf();
    }
    else
    {
        phiI = -muI*fvc::snGrad(V)*mesh.magSf();
        phiE = muE*fvc::snGrad(V)*mesh.magSf();
    }

    // Ion density and electron energy density
    for (label corr = 0; corr <= nNonOrth; corr++)
    {
        if (fluxScheme == "scharfetterGummel")
        {
            fvScalarMatrix nEqn(scharfetterGummel(phiI, Di, n) == Sn);
            nEqn.solve();

            fvScalarMatrix wEqn(scharfetterGummel(phiE, De, w) == Sw);
            wEqn.solve();
        }
        else
        {
            fvScalarMatrix nEqn
            (
                fvm::div(phiI, n, "div(phi,n)")
              - fvm::laplacian(Di, n, "laplacian(D,n)")
             == Sn
            );
            nEqn.solve();

            fvScalarMatrix wEqn
            (
                fvm::div(phiE, w, "div(phi,w)")
              - fvm::laplacian(De, w, "laplacian(D,w)")
             == Sw
            );
            wEqn.solve();
        }
    }

    // Errors
    const label nCells = returnReduce(mesh.nCells(), sumOp<label>());
    const scalar h = Foam::sqrt(1.0/nCells);

    // Largest cell Peclet number of the drift, |phi|/(D |Sf| deltaCoeff)
    const surfaceScalarField g
    (
        mesh.magSf()*mesh.deltaCoeffs()
    );
    const scalar PeI = gMax(mag(phiI.internalField())/g.internalField());
    const scalar PeE =
        gMax
        (
            mag(phiE.internalField())
           /(fvc::interpolate(De)().internalField()*g.internalField())
        );

    scalar L2V, LinfV, L2n, Linfn, L2w, Linfw;
    errors(V, mms::Vexact, L2V, LinfV);
    errors(n, mms::nExact, L2n, Linfn);
    errors(w, mms::wExact, L2w, Linfw);

    Info<< "cells " << nCells << ", h " << h
        << ", cell Peclet number ion " << PeI << ", electron " << PeE << nl
        << "error V  L2 " << L2V << "  max " << LinfV << nl
        << "error n  L2 " << L2n << "  max " << Linfn << nl
        << "error w  L2 " << L2w << "  max " << Linfw << nl << endl;

    if (Pstream::master())
    {
        OFstream os(runTime.path()/"mmsErrors.dat");

        os  << "# cells h PeI PeE L2V LinfV L2n Linfn L2w Linfw" << nl
            << nCells << ' ' << h << ' ' << PeI << ' ' << PeE << ' '
            << L2V << ' ' << LinfV << ' ' << L2n << ' ' << Linfn << ' '
            << L2w << ' ' << Linfw << nl;
    }

    // Fields for inspection, in the time directory 1
    runTime.setTime(1.0, 1);

    V.write();
    n.write();
    w.write();

    volScalarField nErr("nError", n);
    setCells(nErr, mms::nExact);
    nErr.internalField() = n.internalField() - nErr.internalField();
    nErr.write();

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
