 #include <iostream> 
 #include <sstream>
 #include "fvCFD.H"
 #include "bound.H"
 #include "multivariateScheme.H"
 #include "zeroGradientFvPatchFields.H"
 #include "patchToPatchInterpolation.H"
 #include "myFunctionLib.H" // some custom functions
 #include "boundaryConditions/plasmaPotential/plasmaPotential.H"
 //#include "boundaryConditions/mappingPotential/mappingPotential.H"
 #include "boundaryConditions/driftDiffusionPositiveIonDensity/driftDiffusionPositiveIonDensity.H"
 #include "boundaryConditions/driftDiffusionElectronDensity/driftDiffusionElectronDensity.H"
 #include "boundaryConditions/driftDiffusionElectronEnergy/driftDiffusionElectronEnergy.H"
 #include "boundaryConditions/electronTemperature/electronTemperature.H"
 #include "boundaryConditions/electronThermalVelocity/electronThermalVelocity.H"
//***************************************************************************************
 int main(int argc, char *argv[])
 {
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMeshAndFields/readSpeciesProperties.H"
	#include "createMeshAndFields/readReactions.H"
	#include "createMeshAndFields/createFields.H"
	#include "createMeshAndFields/createDielectricFields.H"
	#include "createMeshAndFields/readTransportAndReactionCoefficients.H"
	Foam::sleep(1);
	#include "solvePoisson.H"
	#include "interpolateFields.H"
	#include "solveTransportEquations.H"
	while (runTime.run())
	{
		Info<<"Simulation time:: " << runTime.timeName()<<endl;
		runTime++;
		#include "solveTransportEquations.H"
		#include "solvePoisson.H"
		#include "interpolateFields.H"
		#include "adaptiveTime.H"
		runTime.write();
	}
 }
//*****************************************************************************************
