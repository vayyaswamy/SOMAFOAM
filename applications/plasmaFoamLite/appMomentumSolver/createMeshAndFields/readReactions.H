IOdictionary reactions
(
	IOobject
	(
		"reactions",
		runTime.time().constant(),
		runTime.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

//************************************************************************************************************
// there are 5 different variables that we create in this readReactions code.H *******************************
//============================================================================================================
const int reactionCount(readScalar(reactions.lookup("totalNumberOfReactions")));
string reactionName[reactionCount];
string reactionType[reactionCount];
float reactionRate[reactionCount];
float reactionEnergy[reactionCount];
float reactantMatrix[reactionCount][numberOfSpecies+1];
float productMatrix[reactionCount][numberOfSpecies+1];
//************************************************************************************************************
// first we make all entities of both the matrices zero so we dont have garbage values ***********************
//============================================================================================================
for (int aal=0;aal<reactionCount;aal++){
	for (int bbl=0;bbl<numberOfSpecies+1;bbl++){
		reactantMatrix[aal][bbl]=0;productMatrix[aal][bbl]=0; 
	}
}
//***********************************************************************************************************
// Reading each reactions one by one ************************************************************************
//===========================================================================================================
for (int iccc = 0; iccc < reactionCount; iccc++) {
	Foam::List<Foam::string> parsedReaction(reactions.lookup("reaction_"+convertIntToString(iccc+1)));
	std::stringstream iss(parsedReaction[1]);
	int totalWords= 0;int cttto=0; string word2;
	reactionName[iccc]=parsedReaction[0];
	reactionType[iccc]=parsedReaction[2];
	reactionRate[iccc]=convertStringToFloat(parsedReaction[3]);
	reactantMatrix[iccc][speciesList.size()]=convertStringToFloat(parsedReaction[4]);
	reactionEnergy[iccc]=convertStringToFloat(parsedReaction[5]);
	//********************************************************************************************************
	//Lets count how many words are there in the chemical reaction********************************************
	//========================================================================================================
	while (iss >> word2) 
		totalWords=totalWords+1;
	string word[totalWords-1];
	//********************************************************************************************************
	// Now lets iterate each words of the reaction and store them in an array ================================
	// we also count the number of reactants and products ====================================================
	std::stringstream iss2(parsedReaction[1]);
	int countReactants=0, countProducts=0, LHSorRHS=0;
	while (iss2 >> word2) {
		if (word2=="="){
			LHSorRHS=1;
		}
		else{
			word[cttto]=word2;
			cttto=cttto+1;
		}
		if (LHSorRHS==0 && word2 != "="){
			countReactants=countReactants+1;
		}
		if (LHSorRHS==1 && word2 != "="){
			countProducts=countProducts+1;
		}
	}
	countReactants=(countReactants+1)/3;countProducts=(countProducts+1)/3;
	//*********************************************************************************************************
	// Now lets create two 2D arrays (Matrices) named Reactant and Product Matrix------------------------------
	//=========================================================================================================
	for (int loopWords=0; loopWords<totalWords-1;loopWords++){
		for (int speciesLoop=0; speciesLoop<speciesList.size();speciesLoop++){
			if (loopWords < (3*countReactants-1) && speciesList[speciesLoop]==word[loopWords]){
				reactantMatrix[iccc][speciesLoop]=convertStringToFloat(word[loopWords-1]);
			}
			if (loopWords >= (3*countReactants-1)	&& speciesList[speciesLoop]==word[loopWords]){
				productMatrix[iccc][speciesLoop]=convertStringToFloat(word[loopWords-1]);
			}
		}
	}
}
Info << "Read reactions................... Done !!"<< endl << endl;

//*************************************************************************
// Lets print the reaction table created ----------------------------------
//-------------------------------------------------------------------------
// the code below has nothing to do with the simulation. I am just printing 
// the arrays we created above. So that the user can verify that the input-
// parameters in the case folder are written correctly. -------------------
//=========================================================================
Info << "Reactant Matrix--------" << endl;
for (int aal=0;aal<reactionCount;aal++){
	Info << "reaction_" <<aal <<" energy :" <<reactionEnergy[aal] << endl;
	for (int bbl=0;bbl<numberOfSpecies+1;bbl++){
		Info << reactantMatrix[aal][bbl]<< "   " ; 		
	}
	Info << endl ;
}
//=========================================================================
Info << "Product Matrix--------" << endl;
for (int aal=0;aal<reactionCount;aal++){
	for (int bbl=0;bbl<numberOfSpecies+1;bbl++){
		Info << productMatrix[aal][bbl]<< "   " ; 		
	}
	Info << endl ;
}
//=========================================================================