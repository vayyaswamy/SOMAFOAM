 #include <iostream> 
 #include <sstream>
 #include "fvCFD.H"
 #include "bound.H"
 #include "multivariateScheme.H"
 #include "zeroGradientFvPatchFields.H"
 #include "patchToPatchInterpolation.H"
 #include "myFunctionLib.H" 
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
	while (runTime.run())
	{
		Info<<"Simulation time:: " << runTime.timeName()<<endl;
		runTime++;
		#include "interpolateFields.H"
		#include "solvePoisson.H"
		#include "transport/source.H"
		#include "transport/velocity.H"
		#include "transport/particleContinuity.H"
		#include "transport/energyContinuity.H"
		//#include "adaptiveTime.H"
		runTime.write();
	}
 }
//*****************************************************************************************
