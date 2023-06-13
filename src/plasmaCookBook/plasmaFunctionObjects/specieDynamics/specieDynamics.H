/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    specieDynamics

Description

SourceFiles
    specieDynamics.C

\*---------------------------------------------------------------------------*/

#ifndef specieDynamics_H
#define specieDynamics_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "multiSpeciesPlasmaModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class specieDynamics Declaration
\*---------------------------------------------------------------------------*/

class specieDynamics
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- object to write
		wordList objectNames_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specieDynamics
        (
            const specieDynamics&
        );

        //- Disallow default bitwise assignment
        void operator=(const specieDynamics&);


public:

    //- Runtime type information
    TypeName("specieDynamics");


    // Constructors

        //- Construct from components
        specieDynamics
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
