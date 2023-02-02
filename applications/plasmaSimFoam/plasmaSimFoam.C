/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Application
    plasmaSimFoam

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"

	#include "createPlasmaMesh.H"
	#include "createFields.H"

    lduMatrix::debug = 0;

	coupledLduMatrix::debug = 0;

	blockLduMatrix::debug = 0;

	if (solutionDomain == "plasmaDielectric")
	{
		#include "createDielectricMesh.H"
		#include "createDielectricFields.H"

		while (runTime.run())
		{
			#include "detachPatches.H"
			#include "plasmaEqn.H"
			#include "surfaceCharge.H"

		    runTime++;

		    Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
		        << runTime.elapsedCpuTime() << "s" << endl;

			#include "attachPatches.H"
			#include "solvePoissonD.H"

		    if (runTime.write() && restartCapabale)
		    {    
				thermo.Te().write();

			    thermo.T().write();

			    thermo.Tion().write();

			    thermo.p().write();

			    Phi.write();

				E.write();

				forAll(dielectricRegions, i)
				{
					PhiD[i].write();
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
		while (runTime.run())
		{
			// solve plasma equations 
			#include "plasmaEqn.H"

			// postFix increment. Will update after loop is completed 
		    runTime++;

		    Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
		        << runTime.elapsedCpuTime() << "s" << endl;

			// solve poisson equation  
			#include "solvePoisson.H"

	        // set up new time step by matching electron flux of the system to a user specified courant number  
		    scalar Cofactor = mspm().divFe();
		    scalar deltaTNew = MaxCo/(Cofactor+1e-10);
		    deltaTNew = min(deltaTNew,deltaTMax);
		    deltaTNew = max(deltaTNew,deltaTMin);
		    runTime.setDeltaT(deltaTNew);

		    // print obtained time step and courant number into terminal 
		    Info << "New timestep = " << runTime.deltaTValue() << endl;
		    Info << "Courant = " << Cofactor*runTime.deltaTValue() << endl;

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
						mspm().N(i),
						Y[i].boundaryField().types()
					);
					specN.write();
				}
		    }
		}
	}

    return(0);
}


// ************************************************************************* //
