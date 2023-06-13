/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*----------------------------------------------------------------------------*/

#include "electricalCharacteristics.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electricalCharacteristics, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        electricalCharacteristics,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electricalCharacteristics::electricalCharacteristics
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
    probeFilePtr_(NULL),
	writeFrequency_(-1),
	cellCount_(-1)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

	dict.lookup("writeFrequency") >> writeFrequency_;

	cellCount_ = mesh.cells().size();

    if (mesh.nSolutionD() != 1)
    {
        FatalErrorIn("electricalCharacteristics::")
            << "only 1D meshes are supported.  "
                << abort(FatalError);
    }

    if (Pstream::master())
    {
        if (!time_.processorCase())
        {
            mkDir
            (
                time_.path()
               /"spatioTemporalProbe"
            );

            probeFilePtr_ =
                new OFstream
                (
                    time_.path()
                   /"spatioTemporalProbe"
                   /"spatioTemporalProbe_"+time_.timeName()+".dat"
                );
        }
        else
        {
            mkDir
            (
                time_.path()/".."/"spatioTemporalProbe"
            );

            probeFilePtr_ =
                new OFstream
                (
                    time_.path()/".."
                   /"spatioTemporalProbe"
                   /"spatioTemporalProbe_"+time_.timeName()+".dat"
                );
        }

        (*probeFilePtr_)
            << "VARIABLES=\"t\" \"x\" \"ePower\" \n";
		(*probeFilePtr_)
			 << "ZONE I="<< cellCount_ << ", J=" << writeFrequency_ << ", F=POINT \n";
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::electricalCharacteristics::start()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

	const multiSpeciesPlasmaModel& mspm
	    = mesh.lookupObject<multiSpeciesPlasmaModel>("plasmaProperties");

	const vectorField& xCoord = mesh.C().internalField();

	word Eref = "E";

	word electron = "electron";

	const volVectorField& eField = mesh.lookupObject<volVectorField>
	(
	    Eref
	);

    volScalarField ePowerDeposition
    (
        IOobject
        (
            "ePowerDeposition",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        plasmaConstants::eCharge*(mspm.J(electron) & eField)
    );

    if (Pstream::master())
    {
		if (mesh.time().outputTime())
		{
		    if (!time_.processorCase())
		    {
		        probeFilePtr_ =
		            new OFstream
		            (
		                time_.path()
		               /"spatioTemporalProbe"
		               /"spatioTemporalProbe_"+time_.timeName()+".dat"
		            );
		    }
		    else
		    {
		        probeFilePtr_ =
		            new OFstream
		            (
		                time_.path()/".."
		               /"spatioTemporalProbe"
		               /"spatioTemporalProbe_"+time_.timeName()+".dat"
		            );
		    }
	
			(*probeFilePtr_)
				 << "ZONE I="<< cellCount_ << ", J=" << writeFrequency_ << ", F=POINT \n";
		}

        probeFilePtr_->precision(14);

		forAll(ePowerDeposition, cellI)
		{
		    (*probeFilePtr_) << time_.value() << tab
		        << xCoord[cellI].component(0) << tab
		        << ePowerDeposition[cellI];
			(*probeFilePtr_) << "\r\n";
		}
        return true;
    }
    return true;
}


bool Foam::electricalCharacteristics::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

	const multiSpeciesPlasmaModel& mspm
	    = mesh.lookupObject<multiSpeciesPlasmaModel>("plasmaProperties");

	const vectorField& xCoord = mesh.C().internalField();

	word Eref = "E";

	word electron = "electron";

	const volVectorField& eField = mesh.lookupObject<volVectorField>
	(
	    Eref
	);

    volScalarField ePowerDeposition
    (
        IOobject
        (
            "ePowerDeposition",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        plasmaConstants::eCharge*(mspm.J(electron) & eField)
    );

	int curTimeIndex_ = time_.timeIndex();

    if (Pstream::master() && (curTimeIndex_ %  writeFrequency_ == 0))
    {
		if (mesh.time().outputTime())
		{
		    if (!time_.processorCase())
		    {
		        probeFilePtr_ =
		            new OFstream
		            (
		                time_.path()
		               /"spatioTemporalProbe"
		               /"spatioTemporalProbe_"+time_.timeName()+".dat"
		            );
		    }
		    else
		    {
		        probeFilePtr_ =
		            new OFstream
		            (
		                time_.path()/".."
		               /"spatioTemporalProbe"
		               /"spatioTemporalProbe_"+time_.timeName()+".dat"
		            );
		    }

		    (*probeFilePtr_)
		        << "VARIABLES=\"t\" \"x\" \"ePower\" \n";
			(*probeFilePtr_)
				 << "ZONE I="<< cellCount_ << ", J=" 
				<< writeFrequency_ << ", F=POINT \n";
		}

        probeFilePtr_->precision(14);

		forAll(ePowerDeposition, cellI)
		{
		    (*probeFilePtr_) << time_.value() << tab
		        << xCoord[cellI].component(0) << tab
		        << ePowerDeposition[cellI];
			(*probeFilePtr_) << "\r\n";
		}
        return true;
    }
    return true;
}


bool Foam::electricalCharacteristics::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

	dict.lookup("writeFrequency") >> writeFrequency_;

    return true;
}

// ************************************************************************* //
