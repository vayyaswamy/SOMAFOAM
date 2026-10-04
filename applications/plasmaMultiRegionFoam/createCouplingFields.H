    // Permittivity seen by the Poisson equation in each region, with the
    // cell value on the boundaries. In the plasma it includes the
    // semi-implicit term; iterativeCoupledPotential reads it by this name
    volScalarField epsEff
    (
        IOobject
        (
            "epsilonEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", eps.dimensions(), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    epsEff.internalField() = eps.internalField();
    epsEff.correctBoundaryConditions();

    PtrList<volScalarField> epsEffD(dielectricRegions.size());

    forAll(dielectricRegions, i)
    {
        epsEffD.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "epsilonEff",
                    runTime.timeName(),
                    dielectricRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dielectricRegions[i],
                dimensionedScalar("zero", epsD[i].dimensions(), 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        epsEffD[i].internalField() = epsD[i].internalField();
        epsEffD[i].correctBoundaryConditions();
    }

    // Interface patches of the plasma on which the potential is fixed and
    // relaxed towards the value the dielectric obtains
    DynamicList<label> couplingPatches;

    forAll(Phi.boundaryField(), patchI)
    {
        if
        (
            isA<iterativeCoupledPotentialFvPatchScalarField>
            (
                Phi.boundaryField()[patchI]
            )
        )
        {
            if
            (
               !refCast<const iterativeCoupledPotentialFvPatchScalarField>
                (
                    Phi.boundaryField()[patchI]
                ).fixesValue()
            )
            {
                FatalErrorIn(args.executable())
                    << "Patch " << mesh.boundary()[patchI].name()
                    << " of the plasma must fix the interface potential"
                    << " (fixesValue yes), and its shadow patch in the"
                    << " dielectric must not" << exit(FatalError);
            }

            couplingPatches.append(patchI);
        }
    }

    if (dielectricRegions.size() && couplingPatches.empty())
    {
        FatalErrorIn(args.executable())
            << "No patch of the plasma potential Phi is of type "
            << iterativeCoupledPotentialFvPatchScalarField::typeName << nl
            << "plasmaMultiRegionFoam solves the regions in turn and needs"
            << " this type on both sides of every plasma-dielectric"
            << " interface (somaFoam uses coupledPotential)."
            << exit(FatalError);
    }

    // Controls of the iteration between the regions, optional dictionary
    // plasmaDielectricCoupling in system/fvSolution
    const dictionary couplingDict
    (
        mesh.solutionDict().subOrEmptyDict("plasmaDielectricCoupling")
    );

    const label couplingMaxIterations
    (
        couplingDict.lookupOrDefault<label>("maxIterations", 50)
    );

    const scalar couplingTolerance
    (
        couplingDict.lookupOrDefault<scalar>("tolerance", 1e-8)
    );

    // Relaxation factor of the interface potential, adapted by Aitken's
    // method and carried from one time step to the next
    scalar couplingRelaxation
    (
        couplingDict.lookupOrDefault<scalar>("initialRelaxation", 0.5)
    );

    Info<< "Plasma-dielectric iteration: " << couplingPatches.size()
        << " interface patches, tolerance " << couplingTolerance
        << ", at most " << couplingMaxIterations << " iterations" << endl;
