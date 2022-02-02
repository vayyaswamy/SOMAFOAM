/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

InClass
    Foam::makeMultiSpeciesPlasmaModel

Description
    Macros for instantiating multispeciesplasmamodel

\*---------------------------------------------------------------------------*/


#ifndef makeMultiSpeciesPlasmaModel_H
#define makeMultiSpeciesPlasmaModel_H

#include "addToRunTimeSelectionTable.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

#define makeMultiSpeciesPlasmaModel(PlasmaModel,Thermo)                 		\
                                                                              	\
typedef PlasmaModel<Thermo>                                                		\
    PlasmaModel##Thermo;                                                   		\
                                                                              	\
defineTemplateTypeNameAndDebugWithName                                        	\
(                                                                             	\
    PlasmaModel##Thermo,                                                   		\
    #PlasmaModel                                                           		\
        "<"#Thermo">",                                                        	\
    0                                                                         	\
);                                                                            	\
                                                                              	\
addToRunTimeSelectionTable                                                    	\
(                                                                             	\
    multiSpeciesPlasmaModel,                                               		\
    PlasmaModel##Thermo,                                                   		\
    fvMesh                                                                    	\
)

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
