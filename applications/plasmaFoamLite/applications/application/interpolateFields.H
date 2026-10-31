EFieldTd= mag(EField)/backgroundGasDensity/dimensionCorrection/1E-21; //eField in townsends unit
eTemp=eDensity/(3/2)/kBoltz/(nDensity[0]+cosmicRadiation*dimensionCorrection); // its statistical phys. particles <> 0
//boundMinMax(eTemp, 300*dimensionCorrectionTemp, 220000*dimensionCorrectionTemp); //upper and lower limit of eTemperature
//eDensity=(3/2)*kBoltz*(nDensity[0]+1*dimensionCorrection)*eTemp; //
//-------------------------------------------------------------------------------------------------------------------
// ---------------------------------- mobility interpolation --------------------------------------------------------
forAll (speciesList, intp){
    if (speciesProperties[intp][2]== "EON")
        interpolation(EFieldTd, mobility[intp], mobilityTable,intp,mobRowSizeList[intp] );
    if (speciesProperties[intp][2]== "eTemp")
        interpolation(eTemp, mobility[intp], mobilityTable,intp,mobRowSizeList[intp] ); 
    if (speciesProperties[intp][2]== "einstien"){
        Info << "we dont currently allow einstien relation on mobility" << endl ;
        exit(0);
    }      
    mobility[intp].correctBoundaryConditions();
}
//-------------------------------------------------------------------------------------------------------------------
// ---------------------------------- diffusion interpolation -------------------------------------------------------
forAll (speciesList, intp){
    if (speciesProperties[intp][3] == "EON")
        interpolation(EFieldTd, diffusion[intp], diffusionTable,intp,diffRowSizeList[intp] );
    if (speciesProperties[intp][3]== "eTemp")
        interpolation(eTemp, diffusion[intp], diffusionTable,intp,diffRowSizeList[intp] );
    if (speciesProperties[intp][3]== "einstien"){
        if (intp==0)//electron is always particle 0 in this software
            diffusion[intp]=mobility[intp]*kBoltz*eTemp/(eCharge);
        else //for ions
            diffusion[intp]=mobility[intp]*kBoltz*gasTemperature*dimensionCorrectionTemp/(eCharge);
    }
    diffusion[intp].correctBoundaryConditions();
}
//-------------------------------------------------------------------------------------------------------------------
// ---------------------------------- reaction interpolation --------------------------------------------------------
for(int intp=0;intp<reactionCount ; intp++){
    if (reactionType[intp] == "EON")
        interpolation(EFieldTd, reaction[intp], reactionTable,intp,reactionRowSizeList[intp] );
    if (reactionType[intp]== "eTemp")
        interpolation(eTemp, reaction[intp], reactionTable,intp,reactionRowSizeList[intp] );
    reaction[intp].correctBoundaryConditions();
}
