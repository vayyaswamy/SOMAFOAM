/*---------------------------------------------------------------------------*\
License
    This file is part of SOMAFOAM, distributed under the GNU General Public
    License, version 3 or later.

Application
    icpFieldTest

Description
    Verification of the field solve of the inductive (ICP) option of
    somaFoam (examples/verification/icpSkinEffect). It builds the coil current
    and solves for the azimuthal vector potential with the functions of
    applications/somaFoam/inductiveField.H, which somaFoam uses, on the mesh
    of the case, with the cold-plasma conductivity of an analytic electron
    density profile instead of a plasma solution. system/icpFieldDict:

        frequency   13.56e6;        // [Hz]
        n0          1e17;           // [1/m3]
        nuM         1e7;            // [1/s], momentum-transfer collisions
        profile     uniform;        // n0 for r < R (uniform) or
                                    // n0 J0(2.405 r/R) sin(pi z/L) (besselSine)
        R           0.05;           // [m]
        L           0.05;           // [m] (besselSine)
        coilW       2.5e-3;         // [m], side of the square conductors
        coilTurns   ((0.07 0.00125) (0.07 0.00375) ...);    // (r z) [m]

    The fields AR and AI for 1 A per turn are written to the time 1, and the
    power per square ampere to the log.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "inductiveField.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "icpFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar f = readScalar(dict.lookup("frequency"));
    const scalar n0 = readScalar(dict.lookup("n0"));
    const scalar nuM = readScalar(dict.lookup("nuM"));
    const word profile(dict.lookupOrDefault<word>("profile", "uniform"));
    const scalar R = readScalar(dict.lookup("R"));
    const scalar L = dict.lookupOrDefault<scalar>("L", GREAT);
    const scalar coilW = readScalar(dict.lookup("coilW"));
    const List<vector2> turns(inductiveField::coilTurns(dict));

    const scalar w = 2*mathematicalConstant::pi*f;
    const scalar mu0 = 4e-7*mathematicalConstant::pi;
    const scalar qe = 1.602176634e-19, me = 9.1093837e-31;

    volScalarField AR
    (
        IOobject("AR", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField AI
    (
        IOobject("AI", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );

    // Block field from the boundary types of AR, as in createInductiveFields.H
    wordList aTypes(AR.boundaryField().size());
    forAll(aTypes, p)
    {
        const fvPatchScalarField& pf = AR.boundaryField()[p];
        const word& t = pf.type();
        aTypes[p] =
            (t == "wedge" || t == "empty" || t == "processor" || t == "symmetryPlane")
          ? t
          : (pf.fixesValue() ? word("fixedValue") : word("zeroGradient"));
    }
    volVector2Field A
    (
        IOobject("A", runTime.timeName(), mesh),
        mesh,
        dimensioned<vector2>("zero", AR.dimensions(), vector2::zero),
        aTypes
    );

    // Scaled conductivity, as in solveInductive.H
    volScalarField sigR
    (
        IOobject("sigR", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimless/dimArea, 0.0)
    );
    volScalarField sigI(sigR);

    const vectorField& C = mesh.C().internalField();
    forAll(C, k)
    {
        const scalar r = Foam::sqrt(sqr(C[k].x()) + sqr(C[k].y()));
        const scalar z = C[k].z();

        scalar n = 0;
        if (profile == "uniform")
        {
            n = (r < R) ? n0 : 0;
        }
        else if (r < R && z < L)
        {
            const scalar x = 2.405*r/R;
            scalar j0 = 0, term = 1;
            for (label m = 0; m < 20; m++)
            {
                j0 += term;
                term *= -sqr(x/2)/sqr(scalar(m + 1));
            }
            n = max(n0*j0*Foam::sin(mathematicalConstant::pi*z/L), 0.0);
        }

        const scalar s0 = n*qe*qe/me/(sqr(nuM) + sqr(w));
        sigR[k] = w*mu0*s0*nuM;
        sigI[k] = -w*mu0*s0*w;
    }

    scalarField J;
    const scalar Itotal =
        inductiveField::coilCurrentDensity(mesh, turns, coilW, J);

    Info<< "cells " << returnReduce(mesh.nCells(), sumOp<label>())
        << ", coil current " << Itotal << " A for " << turns.size()
        << " turns" << endl;

    BlockSolverPerformance<vector2> perf =
        inductiveField::solve(mesh, sigR, sigI, J, AR, AI, A);

    // Power per square ampere, 1/2 Re(sigma) w^2 |A|^2, over the wedge
    const scalarField& V = mesh.V();
    scalar P1 = 0;
    forAll(C, k)
    {
        P1 += 0.5*sigR[k]/(w*mu0)*sqr(w)
            *(sqr(AR[k]) + sqr(AI[k]))*V[k];
    }
    reduce(P1, sumOp<scalar>());

    Info<< "block solve iterations " << perf.nIterations()
        << ", power per (1 A)^2 in the wedge " << P1 << " W" << endl;

    runTime.setTime(1.0, 1);
    AR.write();
    AI.write();

    Info<< "End" << endl;

    return 0;
}
