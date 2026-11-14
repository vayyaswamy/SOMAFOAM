// Converged outer loop with step rejection for solutionDomain plasma, as in
// the plasmaDielectric branch of somaFoam.C (see the description there):
// outer iterations (plasmaOuterIteration.H) until the charged species stop
// changing (chargedSpeciesTolerance), then one final iteration; optionally
// the electron energy equation in every iteration
// (electronTemperatureInOuterLoop), whose energy density Ne*Te must then
// also stop changing. A step that does not converge within
// nOuterCorrectors is undone (all registered fields, the surface charge of
// the dielectric side walls included) and retried with half the step,
// down to deltaTMin; the step after a rejection grows by at most
// deltaTIncreaseAfterRejection.
//
// With electronTemperatureInOuterLoop, TeInnerTolerance (PIMPLE, default 0
// = off) also requires the inner passes of the electron energy equation to
// end with a relative change below it: an iteration whose Te passes do not
// converge does not count as converged, and a change above 1 (or not a
// number) ends the iterations at once, so the step is retried before the
// electron energy equation diverges.
//
// Needs createStepControls.H, pimpleDict and chargedSpeciesTolerance, and
// sets TeSolvedInOuterLoop and CofactorTe (step factor of the last
// electron energy solve).

{
    const label nOuterMax =
        pimpleDict.lookupOrDefault<label>("nOuterCorrectors", 1);

    const Switch stepRejection
    (
        pimpleDict.lookupOrDefault<Switch>("stepRejection", true)
    );

    const Switch electronTemperatureInOuterLoop
    (
        pimpleDict.lookupOrDefault<Switch>
        (
            "electronTemperatureInOuterLoop",
            false
        )
    );

    const scalar TeInnerTolerance =
        pimpleDict.lookupOrDefault<scalar>("TeInnerTolerance", 0);

    const bool solveTe = electronTemperatureInOuterLoop && !mspm().fixesTe();

    const label electronI = mspm().species()["electron"];

    TeSolvedInOuterLoop = electronTemperatureInOuterLoop;

    const scalar tStart = runTime.value() - runTime.deltaTValue();

    autoPtr<stepSnapshot> snapshot;

    if (stepRejection)
    {
        snapshot.reset(new stepSnapshot(mesh));
        snapshot().store();
    }

    bool stepRejected = false;
    label outerI = 0;
    scalar lastOuterChange = 0;
    scalar lastTeInner = -1;

    while (true)
    {
        outerI = 0;
        bool finalOuter = false;
        bool outerConverged = false;
        bool outerDiverged = false;

        while (true)
        {
            outerI++;

            if (finalOuter || outerI >= nOuterMax)
            {
                finalOuter = true;
                mesh.solutionDict().add("finalIteration", true);
            }

            #include "plasmaOuterIteration.H"

            scalar change = mspm().chargedSpeciesResidual();

            bool TeInnerConverged = true;
            bool TeInnerDiverged = false;

            if (solveTe)
            {
                const scalarField NeTePrev
                (
                    mspm().N(electronI).internalField()
                   *thermo.Te().internalField()
                );

                {
                    const double tTimerE = solverTimers::now();

                    CofactorTe = pem.ecorrect(chemistry, E);

                    solverTimers::energy += solverTimers::now() - tTimerE;
                }

                const scalarField NeTe
                (
                    mspm().N(electronI).internalField()
                   *thermo.Te().internalField()
                );

                change = max
                (
                    change,
                    gMax(mag(NeTe - NeTePrev))/max(gMax(NeTe), VSMALL)
                );

                if (TeInnerTolerance > 0)
                {
                    lastTeInner = mspm().lastInnerChange("Te");

                    TeInnerConverged = lastTeInner <= TeInnerTolerance;
                    TeInnerDiverged = lastTeInner > 1 || lastTeInner != lastTeInner;
                }
            }

            lastOuterChange = change;

            if (finalOuter)
            {
                mesh.solutionDict().remove("finalIteration");
                break;
            }

            if
            (
                outerI > 1 && change < chargedSpeciesTolerance
             && TeInnerConverged
            )
            {
                finalOuter = true;
                outerConverged = true;
            }
            else if
            (
                stepRejection
             && (
                    (outerI > 2 && (change > 1 || change != change))
                 || TeInnerDiverged
                )
            )
            {
                outerDiverged = true;
                break;
            }
        }

        if (outerConverged)
        {
            Info<< "Outer iterations: " << outerI << ", converged" << endl;

            if (maxIonizationNumber <= 0)
            {
                break;
            }

            {
                const scalarField nu
                (
                    chemistry.dRRDi(ionizationElectronI)().internalField()
                );
                const scalarField& Ne =
                    mspm().N(ionizationElectronI).internalField();
                const scalarField& TeCells = thermo.Te().internalField();
                const scalar TeLimit = 0.999*mspm().TeMax();

                scalar nuMax = 0;
                forAll(nu, celli)
                {
                    if
                    (
                        Ne[celli] >= ionizationNumberMinDensity
                     && TeCells[celli] < TeLimit
                    )
                    {
                        nuMax = max(nuMax, nu[celli]);
                    }
                }
                nuIonizationMax = returnReduce(nuMax, maxOp<scalar>());
            }

            const scalar deltaT = runTime.deltaTValue();
            const scalar ionNumber = nuIonizationMax*deltaT;

            if
            (
                ionNumber <= maxIonizationNumber
             || !stepRejection
             || deltaT <= deltaTMin*(1 + SMALL)
            )
            {
                break;
            }

            snapshot().restore();

            const scalar deltaTRetry = max
            (
                min(0.5*deltaT, 0.9*maxIonizationNumber/nuIonizationMax),
                deltaTMin
            );

            runTime.setDeltaT(deltaTRetry, false);
            runTime.setTime(tStart + deltaTRetry, runTime.timeIndex());

            stepRejected = true;

            Info<< "Step rejected: ionization number " << ionNumber
                << " > " << maxIonizationNumber
                << ", retrying with deltaT = " << deltaTRetry
                << ", Simulation Time = " << runTime.timeName() << "s" << endl;

            continue;
        }

        Info<< "Outer iterations: " << outerI
            << ", NOT converged: relative change " << lastOuterChange
            << " (tolerance " << chargedSpeciesTolerance << ")";
        if (TeInnerTolerance > 0 && solveTe)
        {
            Info<< ", Te inner change " << lastTeInner
                << " (tolerance " << TeInnerTolerance << ")";
        }
        Info<< endl;

        const scalar deltaT = runTime.deltaTValue();

        if (!stepRejection)
        {
            break;
        }

        if (deltaT <= deltaTMin*(1 + SMALL))
        {
            if (outerDiverged)
            {
                FatalErrorIn(args.executable())
                    << "Outer iterations diverge at deltaTMin = "
                    << deltaTMin << ", time " << runTime.value()
                    << exit(FatalError);
            }

            Info<< "Step accepted unconverged at deltaTMin" << endl;
            break;
        }

        snapshot().restore();

        const scalar deltaTRetry = max(0.5*deltaT, deltaTMin);

        runTime.setDeltaT(deltaTRetry, false);
        runTime.setTime(tStart + deltaTRetry, runTime.timeIndex());

        stepRejected = true;

        Info<< "Step rejected: retrying with deltaT = " << deltaTRetry
            << ", Simulation Time = " << runTime.timeName() << "s" << endl;
    }

    if (stepRejected)
    {
        deltaTRejectionCap = runTime.deltaTValue();
    }
    else if (deltaTRejectionCap < GREAT && 2*outerI <= nOuterMax)
    {
        deltaTRejectionCap *= pimpleDict.lookupOrDefault<scalar>
        (
            "deltaTIncreaseAfterRejection",
            1.2
        );

        if (deltaTRejectionCap >= deltaTMax)
        {
            deltaTRejectionCap = GREAT;
        }
    }
}
