/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    electricalCharacteristics

Description

SourceFiles
    electricalCharacteristics.C

\*---------------------------------------------------------------------------*/

#ifndef electricalCharacteristics_H
#define electricalCharacteristics_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "multiSpeciesPlasmaModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class electricalCharacteristics Declaration
\*---------------------------------------------------------------------------*/

class electricalCharacteristics
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

        //- Output file stream
        OFstream* probeFilePtr_;

		int writeFrequency_;

		int cellCount_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        electricalCharacteristics
        (
            const electricalCharacteristics&
        );

        //- Disallow default bitwise assignment
        void operator=(const electricalCharacteristics&);

public:

    //- Runtime type information
    TypeName("electricalCharacteristics");


    // Constructors

        //- Construct from components
        electricalCharacteristics
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
