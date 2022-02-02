/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "emPlasmaFields.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(emPlasmaFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emPlasmaFields::emPlasmaFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    objectNames_()
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "emPlasmaFields::emPlasmaFields"
            "(const objectRegistry&, const dictionary&)"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::emPlasmaFields::~emPlasmaFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::emPlasmaFields::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("objectNames") >> objectNames_;
    }
}


void Foam::emPlasmaFields::execute()
{
    // Do nothing - only valid on write
}


void Foam::emPlasmaFields::end()
{
    // Do nothing - only valid on write
}


void Foam::emPlasmaFields::write()
{
    if (active_)
    {
        forAll(objectNames_, i)
        {
		    if (objectNames_[i] == "ePowerDeposition")
		    {
				word Eref = "E";

				word electron = "electron";

				const volVectorField& E = obr_.lookupObject<volVectorField>
				(
				    Eref
				);

				const multiSpeciesPlasmaModel& mspm
				    = obr_.lookupObject<multiSpeciesPlasmaModel>("plasmaProperties");

		        Info<< "Calculating electron power deposition field." << endl;

		        volScalarField ePowerDeposition
		        (
		            IOobject
		            (
		                "ePowerDeposition",
		                obr_.time().timeName(),
		                obr_,
		                IOobject::NO_READ
		            ),
		            plasmaConstants::eCharge*(mspm.J(electron) & E)
		        );

		        ePowerDeposition.write();

		        Info<< "emPlasmaFields written." << nl << endl;
		    }
		    else
		    {
		        Info<< "Object name (" << objectNames_[i]
		            << ") is not present/implemented. "
		            << endl;
		    }
		}
    }
}


// ************************************************************************* //
