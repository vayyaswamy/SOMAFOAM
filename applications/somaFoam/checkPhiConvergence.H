{
    scalarField magPhi = mag(Phi.internalField());

    //Info << "Phiw = " << mag(Phi.boundaryField()) << endl;

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
    //Info << "relativePhiResidual = " << relativePhiResidual << endl;
}

