/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Application
    somaFoam

Description
    plasma/dielectric
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "coupledFvMatrices.H"
#include "regionCouplePolyPatch.H"
#include "multivariateScheme.H"
#include "regionProperties.H"
#include "zeroGradientFvPatchFields.H"
#include "multiSpeciesPlasmaModel.H"
#include "plasmaEnergyModel.H"
#include "thermoPhysicsTypes.H"
#include "emcModels.H"
#include "pimpleControl.H"
#include "dynamicFvMesh.H"
#include "magLongDelta.H"
#include "staticFvMesh.H"
#include "plasmaPotential.H"
#include "currentTimeStepControl.H"
#include "OFstream.H"
#include "stepSnapshot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "runElectronBoltzmann.H"

    #include "createPlasmaMesh.H"


    #include "createFields.H"

    pimpleControl pimple(mesh);


    if (solutionDomain == "plasmaDielectric")
    {
        // The plasma-dielectric coupling needs a fixed mesh, except for 1D
        // refinement of the plasma mesh (dynamicRefine1DFvMesh): it never
        // splits the boundary faces, so the interface to the dielectric
        // stays the same. The dielectric mesh is not refined.
        if (!isA<staticFvMesh>(mesh) && mesh.type() != "dynamicRefine1DFvMesh")
        {
            FatalErrorIn(args.executable())
                << "A dynamic mesh (" << mesh.type() << ") is selected in"
                << " constant/dynamicMeshDict, but mesh changes are only"
                << " supported for solutionDomain plasma (and dynamicRefine1DFvMesh"
                << " for plasmaDielectric)." << nl
                << "Remove constant/dynamicMeshDict or select staticFvMesh."
                << exit(FatalError);
        }


        #include "createDielectricMesh.H"
        #include "createDielectricFields.H"
        #include "createTotalCurrentFields.H"

        // Set when the plasma mesh changes (1D refinement)
        bool dielectricWeightsOutdated = false;

        // Upper limit of the step after a rejected step (see below); GREAT
        // when none
        scalar deltaTRejectionCap = GREAT;

        while (runTime.run())
        {

            runTime++;

            Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
                    << runTime.elapsedCpuTime() << "s" << endl;

            // Adaptive 1D refinement of the plasma mesh, with the plasma side
            // of the interface detached (as for the plasma equations)
            #include "prolongStore.H"

            if (mesh.update())
            {
                #include "prolongApply.H"

                // The interface face is never split, but the plasma cell
                // behind it can be: the dielectric interpolation weights and
                // distances, which use the cell centres across the
                // interface, are recalculated after the next attach (below)
                dielectricWeightsOutdated = true;

                Info<< "Mesh changed: "
                    << returnReduce(mesh.nCells(), sumOp<label>())
                    << " cells" << endl;
            }

            // Outer iterations. With chargedSpeciesTolerance in PIMPLE they
            // continue until the charged species stop changing: the largest
            // change of a charged density in an iteration, divided by the
            // largest total charged density sum |z| N, below the tolerance
            // (multiSpeciesPlasmaModel::chargedSpeciesResidual). One more,
            // final, iteration then follows (relaxation factors
            // "<field>Final"), as in pimpleControl. pimpleControl's own
            // residualControl takes the last solve of each field, which for
            // species with inner passes is always converged, and normalised
            // residuals stay near 1 for species near the density floor, so
            // neither can judge this. Reactions between charged species
            // conserve charge only when all species are converged to one
            // state.
            //
            // Step rejection (stepRejection, default on with
            // chargedSpeciesTolerance): a step whose outer iterations do not
            // converge within nOuterCorrectors is undone (all registered
            // fields back to the start of the step) and retried with half
            // the step, down to deltaTMin. A relative change above 1 after
            // the second iteration (larger than the total charged density)
            // ends the iterations at once. At deltaTMin an unconverged step
            // is accepted with a warning; a diverging one stops the run.
            // The retry stays at the same time index, so the old-time fields
            // keep the start of the step (Euler time stepping only).
            const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

            const scalar chargedSpeciesTolerance =
                pimpleDict.lookupOrDefault<scalar>("chargedSpeciesTolerance", 0);

            bool stepRejected = false;

            // Set when the outer loop has solved the electron energy equation
            bool TeSolvedInOuterLoop = false;

            if (chargedSpeciesTolerance > 0)
            {
                const label nOuterMax =
                    pimpleDict.lookupOrDefault<label>("nOuterCorrectors", 1);

                const Switch stepRejection
                (
                    pimpleDict.lookupOrDefault<Switch>("stepRejection", true)
                );

                // Optional (electronTemperatureInOuterLoop, default off):
                // electron energy equation in every outer iteration, so that
                // the species fluxes (mobility, rates) and the electron
                // temperature converge together. By default it is solved
                // once after the loop and the species use the electron
                // temperature of the start of the step, a splitting error of
                // first order in the step. In the 1 MHz air DBD both converge
                // to the same pulse as the step goes to zero; in the loop the
                // step costs 2-3 times more and is not more accurate at large
                // steps. With it on, the loop also needs the electron energy
                // density Ne*Te to stop changing: largest change in a cell
                // over its largest value (cells without electrons, whose Te
                // is not determined, do not hold it up).
                const Switch electronTemperatureInOuterLoop
                (
                    pimpleDict.lookupOrDefault<Switch>
                    (
                        "electronTemperatureInOuterLoop",
                        false
                    )
                );

                const label electronI = mspm().species()["electron"];

                TeSolvedInOuterLoop = electronTemperatureInOuterLoop;

                const scalar tStart = runTime.value() - runTime.deltaTValue();

                PtrList<stepSnapshot> snapshots(dielectricRegions.size() + 1);

                if (stepRejection)
                {
                    snapshots.set(0, new stepSnapshot(mesh));

                    forAll(dielectricRegions, i)
                    {
                        snapshots.set
                        (
                            i + 1,
                            new stepSnapshot(dielectricRegions[i])
                        );
                    }

                    forAll(snapshots, i)
                    {
                        snapshots[i].store();
                    }
                }

                label outerI = 0;
                scalar lastOuterChange = 0;

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

                        #include "plasmaDielectricOuterIteration.H"

                        scalar change = mspm().chargedSpeciesResidual();

                        if (electronTemperatureInOuterLoop)
                        {
                            const scalarField NeTePrev
                            (
                                mspm().N(electronI).internalField()
                               *thermo.Te().internalField()
                            );

                            pem.ecorrect(chemistry, E);

                            const scalarField NeTe
                            (
                                mspm().N(electronI).internalField()
                               *thermo.Te().internalField()
                            );

                            change = max
                            (
                                change,
                                gMax(mag(NeTe - NeTePrev))
                               /max(gMax(NeTe), VSMALL)
                            );
                        }

                        lastOuterChange = change;

                        if (finalOuter)
                        {
                            mesh.solutionDict().remove("finalIteration");
                            break;
                        }

                        if (outerI > 1 && change < chargedSpeciesTolerance)
                        {
                            finalOuter = true;
                            outerConverged = true;
                        }
                        else if
                        (
                            stepRejection && outerI > 2
                         && (change > 1 || change != change)
                        )
                        {
                            outerDiverged = true;
                            break;
                        }
                    }

                    if (outerConverged)
                    {
                        Info<< "Outer iterations: " << outerI
                            << ", converged" << endl;
                        break;
                    }

                    Info<< "Outer iterations: " << outerI
                        << ", NOT converged: relative change "
                        << lastOuterChange << " (tolerance "
                        << chargedSpeciesTolerance << ")" << endl;

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

                        Info<< "Step accepted unconverged at deltaTMin"
                            << endl;
                        break;
                    }

                    // Undo the step and retry with half of it
                    forAll(snapshots, i)
                    {
                        snapshots[i].restore();
                    }

                    const scalar deltaTRetry = max(0.5*deltaT, deltaTMin);

                    runTime.setDeltaT(deltaTRetry, false);
                    runTime.setTime(tStart + deltaTRetry, runTime.timeIndex());

                    stepRejected = true;

                    Info<< "Step rejected: retrying with deltaT = "
                        << deltaTRetry << ", Simulation Time = "
                        << runTime.timeName() << "s" << endl;
                }

                // After a rejection the step is held at the accepted one and
                // then grows by at most deltaTIncreaseAfterRejection per step
                // (only on steps that converge within half of
                // nOuterCorrectors), instead of returning at once to the
                // current-limited step that failed
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
            else
            {
                while (pimple.loop())
                {
                    #include "plasmaDielectricOuterIteration.H"
                }
            }

            gradTe = mspm().gradTe();


                scalar Cofactor = mspm().divFe();

                if (!TeSolvedInOuterLoop)
                {
                    pem.ecorrect(chemistry, E);
                }

                #include "postStepTotalCurrent.H"

                scalar deltaTNew = MaxCo/(Cofactor+1e-10);

                // Current-based limit from the electrodes (optional)
                {
                    scalar deltaTI = currentLimitedDeltaT(Phi);

                    forAll(PhiD, i)
                    {
                        deltaTI = min(deltaTI, currentLimitedDeltaT(PhiD[i]));
                    }

                    if (deltaTI < GREAT)
                    {
                        Info<< "Current-limited deltaT = " << deltaTI
                            << endl;
                        deltaTNew = min(deltaTNew, deltaTI);
                    }
                }

                deltaTNew = min(deltaTNew,deltaTMax);

                deltaTNew = min(deltaTNew, deltaTRejectionCap);

                deltaTNew = max(deltaTNew,deltaTMin);

                runTime.setDeltaT(deltaTNew);


            if (runTime.write() && restartCapabale)
            {

                thermo.Te().write();

                thermo.T().write();

                thermo.Tion().write();

                thermo.p().write();

                Phi.write();

                eps.write();


                surfC.write();

                forAll(dielectricRegions, i)
                {
                    PhiD[i].write();

                    ED[i].write();

                    epsD[i].write();


                }

                // Copy N(i) with its boundary conditions, so that their settings
                // (e.g. seec) are written and restarts read them back
                forAll(composition.Y(), i)
                {
                    volScalarField specN
                    (
                        IOobject
                        (
                            composition.species()[i],
                            runTime.timeName(),
                            mesh
                        ),
                        mspm().N(i)
                    );
                    specN.write();
                }
            }
        }
    }
    else if (solutionDomain == "plasma")
    {
        #include "refreshProcessorPatches.H"

        while (runTime.run())
        {


            runTime++;

            Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
                << runTime.elapsedCpuTime() << "s" << endl;

            // Adaptive mesh refinement: the fields are mapped onto the new
            // mesh, and the Poisson equation below is solved on it before
            // the plasma equations use the electric field
            #include "prolongStore.H"

            if (mesh.update())
            {
                #include "prolongApply.H"

                Info<< "Mesh changed: "
                    << returnReduce(mesh.nCells(), sumOp<label>())
                    << " cells" << endl;
            }


            while (pimple.loop())
            {


            while (pimple.correct())
            {
                #include "solvePoisson.H"
            }


            #include "plasmaEqn.H"
            #include "surfaceCharge_new.H"


            }

            gradTe = mspm().gradTe();

            scalar Cofactor1 = mspm().divFe();

        scalar Cofactor2 = pem.ecorrect(chemistry, E);

        scalar Cofactor = max(Cofactor1,Cofactor2);

        meshSize = mspm().meshParameter();


        Info << "Mesh Size Reciprocal = " << gMax(meshSize) << endl;

            scalar deltaTNew = MaxCo/(Cofactor+1e-10);

            // Current-based limit from the electrodes (optional)
            {
                const scalar deltaTI = currentLimitedDeltaT(Phi);

                if (deltaTI < GREAT)
                {
                    Info<< "Current-limited deltaT = " << deltaTI
                        << endl;
                    deltaTNew = min(deltaTNew, deltaTI);
                }
            }

            deltaTNew = min(deltaTNew,deltaTMax);

            deltaTNew = max(deltaTNew,deltaTMin);

            runTime.setDeltaT(deltaTNew);

            Info << "New timestep = " << runTime.deltaTValue() << endl;


            if (runTime.write() && restartCapabale)
            {
                thermo.Te().write();

                thermo.T().write();

                thermo.Tion().write();

                thermo.p().write();

                Phi.write();

                forAll(composition.Y(), i)
                {
                    volScalarField specN
                    (
                        IOobject
                        (
                            composition.species()[i],
                            runTime.timeName(),
                            mesh
                        ),
                        mspm().N(i)
                    );
                    specN.write();
                }
            }
        }
    }
    return(0);
}


// ************************************************************************* //
