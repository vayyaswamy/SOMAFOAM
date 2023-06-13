/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*----------------------------------------------------------------------------*/

#include "specieDynamics.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(specieDynamics, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        specieDynamics,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieDynamics::specieDynamics
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
	objectNames_()
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }
        
	dict.lookup("objectNames") >> objectNames_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::specieDynamics::start()
{
    Info<< "Not implemented." << endl;
    return true;
}


bool Foam::specieDynamics::execute()
{
    Info<< "Not implemented." << endl;
    return true;
}


bool Foam::specieDynamics::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

	dict.lookup("objectNames") >> objectNames_;

    return true;
}

// ************************************************************************* //
