// Change of the potential in the last pass of the Poisson loop
// (solvePoissonD.H), relative to the largest potential in the domain.
//
// It was relative to the potential of each cell, which is near zero next to a
// grounded electrode: there a negligible change gave a relative change of
// order one, the loop never met its tolerance and ran its 100 passes.
{
    const scalar dPhiMax =
        gMax(mag(Phi.internalField() - Phi.prevIter().internalField()));

    scalar PhiScale = gMax(mag(Phi.internalField()));

    forAll(PhiD, i)
    {
        PhiScale = max(PhiScale, gMax(mag(PhiD[i].internalField())));
    }

    relativePhiResidual = dPhiMax/max(PhiScale, SMALL);
}
