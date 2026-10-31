//*****************************************************************************************************
//============ ACCOUNTING FOR THE PRODUCTION OF NEW PARTICLES DUE TO CHEMICAL REACTIONS ===============
//*****************************************************************************************************
eSource=0*eSource;
collFreqxN[0] = 0.0*collFreqxN[0];
for (int ii=0; ii<reactionCount ; ii++){
	reactionSource[ii]=reaction[ii];
	for (int jj=0; jj<numberOfSpecies; jj++){
			reactionSource[ii]=reactionSource[ii]*pow(nDensity[jj],reactantMatrix[ii][jj])
										/pow(dimensionCorrection,reactantMatrix[ii][jj]);
		}
	eSource= eSource + reactionEnergy[ii]*dimensionCorrectionEnergy*reactionSource[ii]/dimensionCorrection;
	
	if (reactionName[ii] == "elastic") //rrbased electron collision freq
		collFreqxN[0] = collFreqxN[0] + reactionSource[ii];
}
// =========================== SOURCE TERM FOR EACH SPECIES ============================================
forAll (speciesList, ii){
	pSource[ii]=0*pSource[ii];
	for (int jj=0;jj<reactionCount;jj++){
		pSource[ii]=pSource[ii]+(productMatrix[jj][ii]-reactantMatrix[jj][ii])*reactionSource[jj];
	}
}// =====================================================================================================