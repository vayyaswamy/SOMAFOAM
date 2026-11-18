// Field solve of the inductive (ICP) option, shared by somaFoam
// (createInductiveFields.H, solveInductive.H) and the verification
// application icpFieldTest (applications/verification/icpFieldTest).
//
// The azimuthal vector potential A = AR + i AI of an axisymmetric em mesh
// (axis along z, r = sqrt(x^2 + y^2)) obeys
//
//     (1/r) d/dr(r dA/dr) + d2A/dz2 - A/r^2 - i w mu0 sigma A = -mu0 J,
//
// written with sigR = w mu0 Re(sigma) and sigI = w mu0 Im(sigma) as the
// block-coupled real system
//
//     lap(AR) - AR/r^2 + sigI AR + sigR AI = -mu0 J
//     lap(AI) - AI/r^2 + sigI AI - sigR AR = 0.

#ifndef inductiveField_H
#define inductiveField_H

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "volVectorNFields.H"

namespace Foam
{
namespace inductiveField
{

//- Current density of the coil for 1 A per turn. Each turn is a square
//  conductor of side coilW centred at (r, z), distributed over the cells by
//  the exact overlap of its cross-section with the (r, z) extent of each
//  cell, so that every turn carries 1 A on any mesh. Returns the total
//  current, which must equal the number of turns.
inline scalar coilCurrentDensity
(
    const fvMesh& emMesh,
    const List<vector2>& turns,
    const scalar coilW,
    scalarField& J
)
{
    const vectorField& Cem = emMesh.C().internalField();
    const pointField& pts = emMesh.points();
    const labelListList& cPts = emMesh.cellPoints();

    J.setSize(Cem.size());
    J = 0;

    scalar total = 0;

    forAll(Cem, i)
    {
        scalar r0 = GREAT, r1 = -GREAT, z0 = GREAT, z1 = -GREAT;

        forAll(cPts[i], k)
        {
            const point& p = pts[cPts[i][k]];
            const scalar r = Foam::sqrt(sqr(p.x()) + sqr(p.y()));
            r0 = min(r0, r);
            r1 = max(r1, r);
            z0 = min(z0, p.z());
            z1 = max(z1, p.z());
        }

        forAll(turns, c)
        {
            const scalar rc = turns[c][0];
            const scalar zc = turns[c][1];

            const scalar lr =
                max(0.0, min(r1, rc + coilW/2) - max(r0, rc - coilW/2));
            const scalar lz =
                max(0.0, min(z1, zc + coilW/2) - max(z0, zc - coilW/2));

            if (lr*lz > 0)
            {
                J[i] += lr*lz/sqr(coilW)/((r1 - r0)*(z1 - z0));
                total += lr*lz/sqr(coilW);
            }
        }
    }

    return returnReduce(total, sumOp<scalar>());
}


//- Solve for (AR, AI) with the source mu0 J and the scaled conductivity
//  (sigR, sigI); A is the block field built from the boundary types of AR
inline BlockSolverPerformance<vector2> solve
(
    const fvMesh& emMesh,
    const volScalarField& sigR,
    const volScalarField& sigI,
    const scalarField& J,
    volScalarField& AR,
    volScalarField& AI,
    volVector2Field& A
)
{
    const scalar mu0 = 4e-7*mathematicalConstant::pi;
    const Time& runTime = emMesh.time();

    volScalarField invR2
    (
        IOobject("invR2", runTime.timeName(), emMesh),
        emMesh,
        dimensionedScalar("zero", dimless/dimArea, 0.0)
    );

    const vectorField& Cem = emMesh.C().internalField();

    forAll(Cem, i)
    {
        invR2.internalField()[i] =
            1.0/max(sqr(Cem[i].x()) + sqr(Cem[i].y()), 1e-12);
    }

    volScalarField src
    (
        IOobject("src", runTime.timeName(), emMesh),
        emMesh,
        dimensionedScalar("zero", AR.dimensions()/dimArea, 0.0)
    );
    src.internalField() = mu0*J;

    fvScalarMatrix eqR
    (
        fvm::laplacian(AR) - fvm::Sp(invR2, AR) + fvm::Sp(sigI, AR) + src
    );
    fvScalarMatrix eqI
    (
        fvm::laplacian(AI) - fvm::Sp(invR2, AI) + fvm::Sp(sigI, AI)
    );
    fvScalarMatrix cRI(fvm::Sp(sigR, AI));
    fvScalarMatrix cIR(fvm::Sp(-sigR, AR));

    fvBlockMatrix<vector2> AEqn(A);
    AEqn.insertEquation(0, eqR);
    AEqn.insertEquation(1, eqI);
    AEqn.insertEquationCoupling(0, 1, cRI);
    AEqn.insertEquationCoupling(1, 0, cIR);

    BlockSolverPerformance<vector2> perf = AEqn.solve();

    AEqn.retrieveSolution(0, AR.internalField());
    AEqn.retrieveSolution(1, AI.internalField());
    AR.correctBoundaryConditions();
    AI.correctBoundaryConditions();

    return perf;
}


//- Turns of the coil from the ICP dictionary: coilTurns ((r z) ...), or the
//  planar coil coilR (r ...) at the height coilZ
inline List<vector2> coilTurns(const dictionary& dict)
{
    if (dict.found("coilTurns"))
    {
        return List<vector2>(dict.lookup("coilTurns"));
    }

    const scalarList coilR(dict.lookup("coilR"));
    const scalar coilZ = readScalar(dict.lookup("coilZ"));

    List<vector2> turns(coilR.size());

    forAll(coilR, c)
    {
        turns[c][0] = coilR[c];
        turns[c][1] = coilZ;
    }

    return turns;
}

} // End namespace inductiveField
} // End namespace Foam

#endif
