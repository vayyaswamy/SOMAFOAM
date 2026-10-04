/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Application
    plasmaMultiRegionFoam

Description
    Plasma with dielectric regions, solved region by region: within each
    time step the Poisson equation of the plasma and the Laplace equations
    of the dielectrics are solved in turn and iterated until the potential
    and the electric displacement match on the interfaces.

    somaFoam (solutionDomain plasmaDielectric) solves the same problem with
    the regions coupled in one matrix. The two solvers use the same case;
    only the type of the potential on the interface patches differs:
    iterativeCoupledPotential here, coupledPotential for somaFoam. The
    plasma step itself is the one of somaFoam.

    Controls of the iteration, optional, in system/fvSolution:

        plasmaDielectricCoupling
        {
            maxIterations       50;
            tolerance           1e-8;   // interface potential change
                                        // relative to the largest potential
            initialRelaxation   0.5;    // then adapted (Aitken)
        }

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
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
#include "iterativeCoupledPotentialFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPlasmaMesh.H"
    #include "createFields.H"

    pimpleControl pimple(mesh);

    if (!isA<staticFvMesh>(mesh))
    {
        FatalErrorIn(args.executable())
            << "A dynamic mesh (" << mesh.type() << ") is selected in"
            << " constant/dynamicMeshDict, but the plasma-dielectric"
            << " coupling assumes a fixed mesh." << exit(FatalError);
    }

    #include "createDielectricMeshes.H"
    #include "createDielectricFields.H"
    #include "createCouplingFields.H"

    while (runTime.run())
    {
        runTime++;

        Info<< "Simulation Time = " << runTime.timeName() << "s" << tab
            << "CPU Time = " << runTime.elapsedCpuTime() << "s" << endl;

        while (pimple.loop())
        {
            while (pimple.correct())
            {
                #include "solvePoissonIterative.H"
            }

            #include "plasmaEqn.H"
            #include "surfaceCharge.H"
        }

        gradTe = mspm().gradTe();

        scalar Cofactor = mspm().divFe();

        pem.ecorrect(chemistry, E);

        scalar deltaTNew = MaxCo/(Cofactor + 1e-10);
        deltaTNew = min(deltaTNew, deltaTMax);
        deltaTNew = max(deltaTNew, deltaTMin);

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

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
