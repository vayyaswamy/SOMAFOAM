/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Application
    refinementLevelsPar

Description
    Transfers the refinement data of an adaptively refined mesh (cellLevel,
    pointLevel and meshModifiers in polyMesh) between the processor meshes
    and the complete mesh, which reconstructParMesh and decomposePar do not
    handle. It is used to decompose a refined mesh again, with balanced
    processor loads, before a parallel run is restarted (see rebalancePar).

    -reconstruct : processor meshes -> complete mesh, after reconstructParMesh
    -decompose   : complete mesh -> processor meshes, after decomposePar

    Both use the cellProcAddressing and pointProcAddressing written by those
    utilities.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "labelIOField.H"
#include "labelIOList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read a list from the polyMesh directory of the most recent mesh
template<class IOListType>
autoPtr<IOListType> readMeshList
(
    const Time& db,
    const word& name,
    const word& instance
)
{
    return autoPtr<IOListType>
    (
        new IOListType
        (
            IOobject
            (
                name,
                instance,
                polyMesh::meshSubDir,
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );
}


// Write a level list into the polyMesh directory of the given instance
void writeLevel
(
    Time& db,
    const word& name,
    const word& instance,
    const labelList& level
)
{
    // Objects are written into the current time directory unless their
    // instance is constant: move the database to the instance of the mesh
    const instant oldTime(db.value(), db.timeName());
    const label oldIndex = db.timeIndex();

    if (instance != db.timeName() && instance != db.constant())
    {
        db.setTime
        (
            instant(readScalar(IStringStream(instance)()), instance),
            oldIndex
        );
    }

    labelIOField levelIO
    (
        IOobject
        (
            name,
            instance,
            polyMesh::meshSubDir,
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelField(level)
    );

    levelIO.write();

    db.setTime(oldTime, oldIndex);
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions(true, false);
    argList::validOptions.insert("reconstruct", "");
    argList::validOptions.insert("decompose", "");

#   include "setRootCase.H"
#   include "createTime.H"

    const bool reconstruct = args.optionFound("reconstruct");
    const bool decompose = args.optionFound("decompose");

    if (reconstruct == decompose)
    {
        FatalErrorIn(args.executable())
            << "Specify either -reconstruct or -decompose"
            << exit(FatalError);
    }

    // Use the last of the selected times
    instantList times = timeSelector::select0(runTime, args);
    runTime.setTime(times[times.size() - 1], times.size() - 1);

    Info<< "Time = " << runTime.timeName() << nl << endl;

    label nProcs = 0;

    while (isDir(args.path()/(word("processor") + name(nProcs))))
    {
        ++nProcs;
    }

    if (!nProcs)
    {
        FatalErrorIn(args.executable())
            << "No processor* directories found"
            << exit(FatalError);
    }

    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );

        databases[procI].setTime(runTime.timeName(), runTime.timeIndex());
    }

    // Instance of the complete mesh
    const word meshInstance
    (
        runTime.findInstance(polyMesh::meshSubDir, "points")
    );

    const wordList levelNames
    (
        IStringStream("(cellLevel pointLevel)")()
    );

    const wordList addressingNames
    (
        IStringStream("(cellProcAddressing pointProcAddressing)")()
    );

    if (reconstruct)
    {
        forAll (levelNames, i)
        {
            dynamicLabelList level;

            forAll (databases, procI)
            {
                Time& db = databases[procI];

                const word procMeshInstance
                (
                    db.findInstance(polyMesh::meshSubDir, "points")
                );

                autoPtr<labelIOField> procLevel
                (
                    readMeshList<labelIOField>
                    (
                        db,
                        levelNames[i],
                        procMeshInstance
                    )
                );

                autoPtr<labelIOList> addr
                (
                    readMeshList<labelIOList>
                    (
                        db,
                        addressingNames[i],
                        procMeshInstance
                    )
                );

                if (procLevel().size() != addr().size())
                {
                    FatalErrorIn(args.executable())
                        << levelNames[i] << " (" << procLevel().size()
                        << ") and " << addressingNames[i] << " ("
                        << addr().size() << ") of processor " << procI
                        << " differ in size in " << procMeshInstance << "."
                        << nl << "Run reconstructParMesh for this time first."
                        << exit(FatalError);
                }

                forAll (addr(), j)
                {
                    const label globalI = addr()[j];

                    if (globalI >= level.size())
                    {
                        level.setSize(globalI + 1, -1);
                    }

                    level[globalI] = procLevel()[j];
                }
            }

            if (findIndex(level, -1) != -1)
            {
                FatalErrorIn(args.executable())
                    << levelNames[i] << " is not defined for all of the "
                    << level.size() << " entries of the complete mesh"
                    << exit(FatalError);
            }

            Info<< "Reconstructed " << levelNames[i] << ": " << level.size()
                << " entries, maximum level " << max(level)
                << ", written to " << meshInstance << endl;

            writeLevel(runTime, levelNames[i], meshInstance, level);
        }

        // The refinement engine settings are the same on all processors
        const fileName procModifiers
        (
            databases[0].path()
           /databases[0].findInstance(polyMesh::meshSubDir, "points")
           /polyMesh::meshSubDir/"meshModifiers"
        );

        if (isFile(procModifiers))
        {
            cp
            (
                procModifiers,
                runTime.path()/meshInstance/polyMesh::meshSubDir
               /"meshModifiers"
            );
        }
    }
    else
    {
        forAll (levelNames, i)
        {
            autoPtr<labelIOField> level
            (
                readMeshList<labelIOField>(runTime, levelNames[i], meshInstance)
            );

            forAll (databases, procI)
            {
                Time& db = databases[procI];

                const word procMeshInstance
                (
                    db.findInstance(polyMesh::meshSubDir, "points")
                );

                autoPtr<labelIOList> addr
                (
                    readMeshList<labelIOList>
                    (
                        db,
                        addressingNames[i],
                        procMeshInstance
                    )
                );

                if (addr().size() && max(addr()) >= level().size())
                {
                    FatalErrorIn(args.executable())
                        << addressingNames[i] << " of processor " << procI
                        << " refers to entry " << max(addr()) << " but "
                        << levelNames[i] << " in " << meshInstance
                        << " has " << level().size() << " entries." << nl
                        << "Run decomposePar for this time first."
                        << exit(FatalError);
                }

                labelList procLevel(addr().size());

                forAll (addr(), j)
                {
                    procLevel[j] = level()[addr()[j]];
                }

                writeLevel(db, levelNames[i], procMeshInstance, procLevel);

                Info<< "Processor " << procI << ": " << levelNames[i] << " "
                    << procLevel.size() << " entries, written to "
                    << procMeshInstance << endl;
            }
        }

        const fileName modifiers
        (
            runTime.path()/meshInstance/polyMesh::meshSubDir/"meshModifiers"
        );

        if (isFile(modifiers))
        {
            forAll (databases, procI)
            {
                cp
                (
                    modifiers,
                    databases[procI].path()
                   /databases[procI].findInstance
                    (
                        polyMesh::meshSubDir,
                        "points"
                    )
                   /polyMesh::meshSubDir/"meshModifiers"
                );
            }
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
