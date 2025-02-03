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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"

	#include "createPlasmaMesh.H"

	//Info << "createFields " << endl;

	#include "createFields.H"

	pimpleControl pimple(mesh);

	//Info << "Done " << endl;

    //lduMatrix::debug = 0;

	//coupledLduMatrix::debug = 0;

	//blockLduMatrix::debug = 0;
	

	if (solutionDomain == "plasmaDielectric")
	{


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


    	   		#include "solvePoissonD.H"
			
				#include "detachPatches.H"

				#include "plasmaEqn.H"
			
				#include "surfaceCharge.H"

				//#include "dielectricJ.H"
		       	
			}	
			
			

		    scalar Cofactor = mspm().divFe();

		    scalar deltaTNew = MaxCo/(Cofactor+1e-10);

		    deltaTNew = min(deltaTNew,deltaTMax);

		    deltaTNew = max(deltaTNew,deltaTMin);

		    runTime.setDeltaT(deltaTNew);

		    //Info << "New timestep = " << runTime.deltaTValue() << endl;

		    //Info << "Courant = " << Cofactor*runTime.deltaTValue() << endl;

		    
    	


		    if (runTime.write() && restartCapabale)
		    {    

				thermo.Te().write();

			    thermo.T().write();

			    thermo.Tion().write();

			    thermo.p().write();

			    Phi.write();

				eps.write();

				//E.write();

				surfC.write();

				forAll(dielectricRegions, i)
				{
					PhiD[i].write();

					ED[i].write();

					epsD[i].write();


				}

			    forAll(composition.Y(), i)
			    {
					/*volScalarField specN
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
					specN.write();*/
					mspm().N(i).write();
			    }
		    }
		}
	}
	else if (solutionDomain == "plasma")
	{
		while (runTime.run())
		{


			runTime++;

			Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
		        << runTime.elapsedCpuTime() << "s" << endl;

		    while (pimple.loop())
		    {
			
			/*while (pimple.correct())
			{
				#include "solvePoisson.H"	
			}*/

		    while (pimple.correct())
		    {
		   		#include "solvePoisson.H" 	
		    }
		    

			#include "plasmaEqn.H"	
			#include "surfaceCharge_new.H"
			
			
			

		    }

		    //pem.ecorrect(chemistry, E);	

	//Info << "Te corrected " << endl;

	


			gradTe = mspm().gradTe();
			
			

		    scalar Cofactor1 = mspm().divFe();
      
        scalar Cofactor2 = pem.ecorrect(chemistry, E);
            
        scalar Cofactor = max(Cofactor1,Cofactor2);

        scalar meshParameter = mspm().meshParameter();

		    Info << "Cofactor1 = " << Cofactor1 << endl;
            
        Info << "Cofactor2 = " << Cofactor2 << endl;
        
        Info << "Cofactor = " << Cofactor << endl;

        Info << "Mesh Parameter = " << meshParameter << endl;

		    scalar deltaTNew = MaxCo/(Cofactor+1e-10);

		    deltaTNew = min(deltaTNew,deltaTMax);

		    deltaTNew = max(deltaTNew,deltaTMin);

		    runTime.setDeltaT(deltaTNew);

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
						mspm().N(i)
					);
					specN.write();
					//mspm().N(i).write();
				}
		    }
		}
	}
    return(0);
}


// ************************************************************************* //
