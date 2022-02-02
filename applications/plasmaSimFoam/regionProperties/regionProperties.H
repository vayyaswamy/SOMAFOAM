/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    regionProperties

Description
    Simple class to hold region information for coupled region simulations

SourceFiles
    regionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef regionProperties_H
#define regionProperties_H

#include "IOdictionary.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public IOdictionary
{
    // Private data

        List<word> dielectricRegionNames_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionProperties(const regionProperties&);

        //- Disallow default bitwise assignment
        void operator=(const regionProperties&);


public:

    // Constructors

        //- Construct from components
        regionProperties(const Time& runTime);


    // Destructor

        ~regionProperties();


    // Member Functions

        // Access

            const List<word>& dielectricRegionNames() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
