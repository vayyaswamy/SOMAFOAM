{
    scalarField magPhi = mag(Phi.internalField());


    forAll(magPhi, cellI)
    {
        if (magPhi[cellI] < SMALL)
        {
            magPhi[cellI] = SMALL;
        }
    }

    relativePhiResidual =
        gMax
        (
            mag
            (
                Phi.internalField()
              - Phi.prevIter().internalField()
            )
           /magPhi
        );
}
