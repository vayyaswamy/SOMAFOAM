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

        // Set when the plasma mesh changes (1D refinement)
        bool dielectricWeightsOutdated = false;

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

            while (pimple.loop())
            {

                #include "attachPatches.H"

                if (dielectricWeightsOutdated)
                {
                    forAll(dielectricRegions, i)
                    {
                        dielectricRegions[i].surfaceInterpolation::movePoints();

                        // Distance between the cell centres across the
                        // interface, used by the harmonic interpolation of
                        // the permittivity (cached per mesh)
                        magLongDelta::Delete(dielectricRegions[i]);
                    }

                    dielectricWeightsOutdated = false;
                }

                while (pimple.correct())
                {
                        #include "solvePoissonD.H"
                }
                #include "detachPatches.H"

                #include "plasmaEqn.H"

                #include "surfaceCharge.H"


            }

            gradTe = mspm().gradTe();


                scalar Cofactor = mspm().divFe();

                scalar Cofactor2 = pem.ecorrect(chemistry, E);

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
