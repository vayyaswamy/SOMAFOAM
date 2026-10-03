/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "electrodeVoltageCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electrodeVoltageCurrent, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        electrodeVoltageCurrent,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::electrodeVoltageCurrent::makeFiles()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    forAll(patchNames_, i)
    {
        if (mesh.boundaryMesh().findPatchID(patchNames_[i]) < 0)
        {
            FatalErrorIn("electrodeVoltageCurrent::makeFiles()")
                << "Patch " << patchNames_[i] << " not found. Valid patches: "
                << mesh.boundaryMesh().names() << exit(FatalError);
        }
    }

    if (!Pstream::master())
    {
        return;
    }

    // A new directory per start time, so a restart does not overwrite
    const fileName outputDir
    (
        (Pstream::parRun() ? time_.path()/".." : time_.path())
       /name_
       /time_.timeName()
    );

    mkDir(outputDir);

    files_.setSize(patchNames_.size());

    forAll(patchNames_, i)
    {
        files_.set(i, new OFstream(outputDir/(patchNames_[i] + ".dat")));

        files_[i].precision(10);

        files_[i]
            << "# patch " << patchNames_[i] << nl
            << "# currents are positive from the electrode into the plasma"
            << nl
            << "# t [s]" << tab << "V [V]" << tab << "I_conduction [A]" << tab
            << "I_displacement [A]" << tab << "I_total [A]" << tab
            << "j_total [A/m2]" << endl;
    }
}


void Foam::electrodeVoltageCurrent::write()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    if
    (
        !mesh.foundObject<volScalarField>(potentialName_)
     || !mesh.foundObject<volVectorField>(conductionCurrentName_)
     || !mesh.foundObject<volVectorField>(totalCurrentName_)
    )
    {
        WarningIn("electrodeVoltageCurrent::write()")
            << "Fields " << potentialName_ << ", " << conductionCurrentName_
            << " or " << totalCurrentName_ << " not found; nothing written"
            << endl;

        return;
    }

    const volScalarField& Phi =
        mesh.lookupObject<volScalarField>(potentialName_);
    const volVectorField& Jcond =
        mesh.lookupObject<volVectorField>(conductionCurrentName_);
    const volVectorField& Jtot =
        mesh.lookupObject<volVectorField>(totalCurrentName_);

    forAll(patchNames_, i)
    {
        const label patchI = mesh.boundaryMesh().findPatchID(patchNames_[i]);

        const vectorField& Sf = mesh.Sf().boundaryField()[patchI];
        const scalarField& magSf = mesh.magSf().boundaryField()[patchI];

        const scalar area = gSum(magSf);

        const scalar V = gSum(Phi.boundaryField()[patchI]*magSf)/area;

        // Sf points out of the plasma, into the electrode
        const scalar Icond = -gSum(Jcond.boundaryField()[patchI] & Sf);
        const scalar Itot = -gSum(Jtot.boundaryField()[patchI] & Sf);

        if (Pstream::master())
        {
            files_[i]
                << time_.value() << tab << V << tab << Icond << tab
                << Itot - Icond << tab << Itot << tab << Itot/area << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrodeVoltageCurrent::electrodeVoltageCurrent
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
    patchNames_(),
    potentialName_("Phi"),
    conductionCurrentName_("Jnet"),
    totalCurrentName_("Jtot"),
    outputInterval_(1),
    files_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::electrodeVoltageCurrent::start()
{
    makeFiles();

    write();

    return true;
}


bool Foam::electrodeVoltageCurrent::execute()
{
    if (time_.timeIndex() % outputInterval_ == 0)
    {
        write();
    }

    return true;
}


bool Foam::electrodeVoltageCurrent::read(const dictionary& dict)
{
    dict.readIfPresent("region", regionName_);
    dict.lookup("patches") >> patchNames_;
    dict.readIfPresent("potential", potentialName_);
    dict.readIfPresent("conductionCurrent", conductionCurrentName_);
    dict.readIfPresent("totalCurrent", totalCurrentName_);
    dict.readIfPresent("outputInterval", outputInterval_);

    if (outputInterval_ < 1)
    {
        outputInterval_ = 1;
    }

    return true;
}


// ************************************************************************* //
