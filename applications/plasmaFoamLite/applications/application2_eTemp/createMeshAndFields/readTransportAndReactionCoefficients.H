//*********************************************************************************************************************
//******** READING MOBILITY DATA AND STORING IN ARRAYS*****************************************************************
//=====================================================================================================================
Info << "Reading transpor and Reaction Coefficient .......................... !! " << endl;
Info << "Reading mobility data!! ***************"<< endl ;
int mobRowSizeList[speciesList.size()];
Info << "--------------------------------------------------------------------------------------" << endl;
Info << "    Number of Species given:: " << speciesList.size() << endl;
Info << "--------------------------------------------------------------------------------------" << endl;
forAll (speciesList, ii){
 if (speciesProperties[ii][2]!="constant" && speciesProperties[ii][2]!="einstien"){
	Info << "         Mobility of "<< speciesList[ii] << " read from a text file" <<endl; 
	mobRowSizeList[ii]=countRows(speciesList[ii]+"_mobility");
 }
 else {
	 Info << "         Mobility of "<< speciesList[ii] << " is set to constant" <<endl; 
 }
}//=====================================================================================================================
float mobilityTable[speciesList.size()][1000][2];
forAll (speciesList,ii){
	if (speciesProperties[ii][2]!="constant" && speciesProperties[ii][2]!="einstien" ){
		string pathOfFile="constant/"+speciesList[ii]+"_mobility";
		std::ifstream iFile(pathOfFile.c_str());
		int cnt=0;
		Info <<endl<< "         " <<speciesList[ii] << " mobility table ---" << endl;
		while(! iFile.eof()){
			iFile >> mobilityTable[ii][cnt][0];iFile >> mobilityTable[ii][cnt][1];
			mobilityTable[ii][cnt][1]=mobilityTable[ii][cnt][1]/backgroundGasDensity;
			Info << "           " <<mobilityTable[ii][cnt][0] << "  " << mobilityTable[ii][cnt][1]<< endl;
			cnt++;
		}
	}
}
Info << "Finished reading mobility data!!*********************************"<< endl<< endl;

//*********************************************************************************************************************
//------------ READING DIFFUSION DATA AND STORING IN ARRAY ============================================================
//=====================================================================================================================
Info << " Reading diffusion data ............!! "<< endl ;
Info << "------------------------------------------------------------------------------" << endl;
Info << "    Number of Species given:: " << speciesList.size() << endl;
Info << "------------------------------------------------------------------------------" << endl;
int diffRowSizeList[speciesList.size()];
forAll (speciesList, ii){
	if (speciesProperties[ii][3]!="constant" && speciesProperties[ii][3]!="einstien" ){
		Info << "         Diffusion of "<< speciesList[ii] << " read from a text file" <<endl; 
		diffRowSizeList[ii]=countRows(speciesList[ii]+"_diffusion");
	}
	else{
		Info << "         Diffusion of "<< speciesList[ii] << " is set to constant" <<endl; 
	}
}
//**********************************************************************************************************************
float diffusionTable[speciesList.size()][1000][2];
forAll (speciesList,ii){
	if (speciesProperties[ii][3]!="constant" && speciesProperties[ii][3]!="einstien" ){
		string pathOfFile="constant/"+speciesList[ii]+"_diffusion";
		std::ifstream iFile(pathOfFile.c_str());
		int cnt=0;
		Info <<endl<< "         " << speciesList[ii] << " diffusion table ---" << endl;
		while(! iFile.eof()){
			iFile >> diffusionTable[ii][cnt][0];iFile >> diffusionTable[ii][cnt][1];
			diffusionTable[ii][cnt][1]=diffusionTable[ii][cnt][1]/backgroundGasDensity;
			Info<< "           " <<diffusionTable[ii][cnt][0] << "  " << diffusionTable[ii][cnt][1]<< endl;
			cnt++;
		}
	}
}
Info << "Finished reading diffusion data!!**********************************" << endl<< endl;


//***********************************************************************************************************************
//------------ READING REACTION RATE DATA AND STORING IN ARRAYS =========================================================
//***********************************************************************************************************************
Info << "Reading reaction rates!! ***************"<< endl ;
int reactionRowSizeList[reactionCount];
Info << "--------------------------------------------------------------------------------------" << endl;
Info << "    Number of reactions given:: " << reactionCount << endl;
Info << "--------------------------------------------------------------------------------------" << endl;
for (int ii=0; ii<reactionCount; ii++){
	if (reactionType[ii]!="constant"){
		Info << "         Reaction_"<< ii+1 << " rate value obtaned from text file" <<endl; 
		string tmpScalarVar="reaction_"+convertIntToString(ii+1);
		reactionRowSizeList[ii]=countRows(tmpScalarVar);
		Info << "rowCount" << reactionRowSizeList[ii];
	}
	else {
		Info << "         Reaction_"<< ii+1 << "  rate value is set to constant" <<endl; 
	}
}
//************************************************************************************************************************
float reactionTable[reactionCount][1000][2];
for (int ii=0; ii<reactionCount; ii++){
	if (reactionType[ii]!="constant"){
		string pathOfFile="constant/reaction_"+convertIntToString(ii+1);
		std::ifstream iFile(pathOfFile.c_str());
		int cnt=0;
		Info<<endl<< "           " << "Reaction_" << (ii+1) << " Reaction table ---" << endl;
		while(! iFile.eof()){
			iFile >> reactionTable[ii][cnt][0];iFile >> reactionTable[ii][cnt][1];
			Info<< "           " <<reactionTable[ii][cnt][0] << "  " << reactionTable[ii][cnt][1]<< endl;
			cnt++;
		}
	}
}
Info << "Finished reading reaction rates!!************************************************************"<< endl<< endl;
Info <<"Code was updated" <<endl;
//*************************************************************************************************************************