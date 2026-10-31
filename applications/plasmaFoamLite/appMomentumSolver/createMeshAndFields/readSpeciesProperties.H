IOdictionary speciesDict
(
	IOobject
	(
		"speciesProperties",
		runTime.time().constant(),
		runTime.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);
const word emSolver(speciesDict.lookup("emSolver"));
//const int emSolver((speciesDict.lookup("emSolver")));

//**************************************************************************************************
Foam::List<Foam::word> speciesList(speciesDict.lookup("speciesList"));
Info << "List of Species" << endl <<" ============================" << endl << speciesList<< endl;
int numberOfSpecies = speciesList.size(); 
Info << "Number of Species::" << numberOfSpecies << endl; 
const float backgroundGasDensity(readScalar(speciesDict.lookup("backgroundGasDensity")));
Info << "Background gas Density ::" << backgroundGasDensity << endl;
const float gasTemperature(readScalar(speciesDict.lookup("gasTemperature")));
Info << "Gas Temperature ::" << gasTemperature << endl;
const float cosmicRadiation(readScalar(speciesDict.lookup("cosmicRadiation")));
Info << "Cosmic Microwave Radiation ::" << cosmicRadiation << endl;
//**************************************************************************************************
// ----- READING PROPERTIES OF EACH SPECIES ========================================================
//**************************************************************************************************
string speciesProperties[numberOfSpecies][6];
string transportModel[numberOfSpecies];
string collisionType[numberOfSpecies];
float chargeNumber[numberOfSpecies];
float speciesMass[numberOfSpecies];
forAll (speciesList, ii){
	Foam::List<Foam::word> parsedSpeciesProperties(speciesDict.lookup(speciesList[ii]));
	speciesProperties[ii][0]=parsedSpeciesProperties[0];// mass
	speciesMass[ii]=convertStringToFloat(parsedSpeciesProperties[0]);// mass
	speciesProperties[ii][1]=parsedSpeciesProperties[1];// charge
	chargeNumber[ii]=convertStringToInt(parsedSpeciesProperties[1]);// charge
	speciesProperties[ii][2]=parsedSpeciesProperties[2];// mobilityModel
	speciesProperties[ii][3]=parsedSpeciesProperties[3];// diffusionModel
	speciesProperties[ii][4]="0";
	speciesProperties[ii][5]="0";
	transportModel[ii] = parsedSpeciesProperties[6];
	collisionType[ii] = parsedSpeciesProperties[7];
	if (parsedSpeciesProperties[2]== "constant"){
		if (parsedSpeciesProperties.size()<5){
			Info << "value of mobility or diffusion missing. Please Check input file"; exit(0);
		}
		speciesProperties[ii][4]=parsedSpeciesProperties[4];}// diffusionModel
	if (parsedSpeciesProperties[3]== "constant"){
		if (parsedSpeciesProperties.size()<6){
			Info << "value of mobility or diffusion missing. Please Check input file"; exit(0);
		}
		speciesProperties[ii][5]=parsedSpeciesProperties[5];
	}	
}

//****************************************************************************************************
// Printing some of the variables created ------------------------------------------------------------
//====================================================================================================
for (int ii=0; ii<speciesList.size();ii++){
	Info <<"Name of species :: " <<speciesList[ii] << "\t Mass ::"<< speciesMass[ii] << 
	"\t Charge :: "<<chargeNumber[ii]<< "\t Mobility model :: " << speciesProperties[ii][2]  
	<< "\t Diffusion Model :: " << speciesProperties[ii][3] ;
	if (speciesProperties[ii][2]=="constant")
		Info << "\t Mobility Value :: " << speciesProperties[ii][4];
	else
		Info << "\t Mobility Value :: from table" ;
	if (speciesProperties[ii][3]=="constant")
		Info << "\t Diffusion Value :: " << speciesProperties[ii][5];
	else
		Info << "\t Diffusion Value :: from table" ;
Info << endl;	
}//====================================================================================================