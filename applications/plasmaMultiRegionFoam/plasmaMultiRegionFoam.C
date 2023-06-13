/*---------------------------------------------------------------------------*\
Copyright (C) 2019 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Application
    plasmaMultiRegionFoam

Description
    plasma-dielectric

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "zeroGradientFvPatchFields.H"
#include "multivariateScheme.H"
#include "multiSpeciesPlasmaModel.H"
#include "plasmaEnergyModel.H"
#include "thermoPhysicsTypes.H"
#include "emcModels.H"
#include "spatioTemporal.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createPlasmaMeshes.H"
    #include "createDielectricMeshes.H"
    #include "plasmaControls.H"
    #include "createPlasmaFields.H"
    #include "createDielectricFields.H"

    lduMatrix::debug = 0;

	blockLduMatrix::debug = 0;

    while(runTime.loop())
    {
		#include "solvePlasma.H"

		#include "solvePoisson.H"

     	if (courantNumber == "on")
	 	{
		    scalar Cofactor = mspm().divFe();
		    scalar deltaTNew = MaxCo/(Cofactor+1e-10);
		    deltaTNew = min(deltaTNew,deltaTMax);
		    deltaTNew = max(deltaTNew,deltaTMin);
		    runTime.setDeltaT(deltaTNew);

		    // print obtained time step and courant number into terminal 
	    	Info << "New timestep = " << runTime.deltaTValue() << endl;
	    	Info << "Courant = " << Cofactor*runTime.deltaTValue() << endl;
	 	}
	 	else
	 	{

	 	}

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

		if (refSpatioTemporal!="noneD")
		{	
			if (refSpatioTemporal == "oneD")
			{
				electronDiv = mspm().divFelectron();
			}

			spatioTemp -> correct();
		}

	    Info<< "Simulation Time = " << runTime.timeName() << "s" << tab << "CPU Time = "
		        << runTime.elapsedCpuTime() << "s" << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
