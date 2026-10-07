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
#include "staticFvMesh.H"

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
        // The plasma-dielectric coupling assumes a fixed mesh
        if (!isA<staticFvMesh>(mesh))
        {
            FatalErrorIn(args.executable())
                << "A dynamic mesh (" << mesh.type() << ") is selected in"
                << " constant/dynamicMeshDict, but mesh changes are only"
                << " supported for solutionDomain plasma." << nl
                << "Remove constant/dynamicMeshDict or select staticFvMesh."
                << exit(FatalError);
        }


        #include "createDielectricMesh.H"
        #include "createDielectricFields.H"

        while (runTime.run())
        {

            runTime++;

            Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
                    << runTime.elapsedCpuTime() << "s" << endl;

            while (pimple.loop())
            {

                #include "attachPatches.H"

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
                        mspm().N(i),
                        Y[i].boundaryField().types()
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
