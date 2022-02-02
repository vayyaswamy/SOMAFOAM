/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.



Application
    plasmaDictionary

Description
    Utility to change dictionary entries, e.g. can be used to change the patch
    type in the field and polyMesh/boundary files.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "volFields.H"
#include "stringListOps.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool merge(dictionary&, const dictionary&, const bool);

// Add thisEntry to dictionary thisDict.
bool addEntry
(
    dictionary& thisDict,
    entry& thisEntry,
    const entry& mergeEntry,
    const bool literalRE
)
{
    bool changed = false;

    // Recursively merge sub-dictionaries
    // TODO: merge without copying
    if (thisEntry.isDict() && mergeEntry.isDict())
    {
        if
        (
            merge
            (
                const_cast<dictionary&>(thisEntry.dict()),
                mergeEntry.dict(),
                literalRE
            )
        )
        {
            changed = true;
        }
    }
    else
    {
        // Should use in-place modification instead of adding
        thisDict.add(mergeEntry.clone(thisDict).ptr(), true);
        changed = true;
    }

    return changed;
}


// Dictionary merging/editing.
// literalRE:
// - true: behave like dictionary::merge, i.e. add regexps just like
//   any other key.
// - false : interpret wildcard as a rule for items to be matched.
bool merge
(
    dictionary& thisDict,
    const dictionary& mergeDict,
    const bool literalRE
)
{
    bool wildCardInMergeDict = false;

    bool changed = false;

    // Save current (non-wildcard) keys before adding items.
    HashSet<word> thisKeysSet;
    {
        List<keyType> keys = thisDict.keys(false);
        forAll(keys, i)
        {
            thisKeysSet.insert(keys[i]);
        }
    }

    // Pass 1. All literal matches

    forAllConstIter(IDLList<entry>, mergeDict, mergeIter)
    {
        const keyType& key = mergeIter().keyword();

        if (literalRE || !key.isPattern())
        {
            entry* entryPtr = thisDict.lookupEntryPtr
            (
                key,
                false,              // recursive
                false               // patternMatch
            );

            if (entryPtr)
            {

                // Mark thisDict entry as having been match for wildcard
                // handling later on.
                thisKeysSet.erase(entryPtr->keyword());

                if
                (
                    addEntry
                    (
                        thisDict,
                       *entryPtr,
                        mergeIter(),
                        literalRE
                    )
                )
                {
                    changed = true;
                }
            }
            else
            {
                // not found - just add
                thisDict.add(mergeIter().clone(thisDict).ptr());
                changed = true;
            }
        }
    }


    // Pass 2. Wildcard matches (if any) on any non-match keys.

    if (!literalRE && thisKeysSet.size() > 0)
    {
        wordList thisKeys(thisKeysSet.toc());

        forAllConstIter(IDLList<entry>, mergeDict, mergeIter)
        {
            const keyType& key = mergeIter().keyword();

            if (key.isPattern())
            {
                // Find all matching entries in the original thisDict

                if (!wildCardInMergeDict)
                {
                    wildCardInMergeDict = true;
                }

                labelList matches = findStrings(key, thisKeys);

                forAll(matches, i)
                {
                    label matchI = matches[i];

                    entry& thisEntry = const_cast<entry&>
                    (
                        thisDict.lookupEntry(thisKeys[matchI], false, false)
                    );

                    if
                    (
                        addEntry
                        (
                            thisDict,
                            thisEntry,
                            mergeIter(),
                            literalRE
                        )
                    )
                    {
                        changed = true;
                    }
                }
            }
        }
    }

    return changed;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("instance", "instance");
    argList::validOptions.insert("literalRE", "");
    argList::validOptions.insert("U", "");
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"


    bool literalRE = args.optionFound("literalRE");
    bool Uref = args.optionFound("U");

    if (literalRE)
    {
        Info<< "Not interpreting any regular expressions (RE)"
            << " in the plasmaDict." << endl
            << "Instead they are handled as any other entry, i.e. added if"
            << " not present." << endl;
    }


    fileName regionPrefix = "";
    if (regionName != fvMesh::defaultRegion)
    {
        regionPrefix = regionName;
//		#include "createRegionFields.H"	
    }
//	else
//	{
//		#include "createFields.H"	
//	}

    word instance = runTime.timeName();
    if (args.options().found("instance"))
    {
        instance = args.options()["instance"];
    }

    // Get the replacement rules from a dictionary
    IOdictionary dict
    (
        IOobject
        (
            "plasmaDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const dictionary& replaceDicts = dict.subDict("dictionaryReplacement");
    Info<< "Read dictionary " << dict.name()
        << " with replacements for dictionaries "
        << replaceDicts.toc() << endl;


    // Every replacement is a dictionary name and a keyword in this

    forAllConstIter(dictionary, replaceDicts, fieldIter)
    {
        const word& fieldName = fieldIter().keyword();
        Info<< "Replacing entries in dictionary " << fieldName << endl;

        // Handle 'boundary' specially:
        // - is PtrList of dictionaries
        // - is in polyMesh/
        if (fieldName == "boundary")
        {
            Info<< "Special handling of " << fieldName
                << " as polyMesh/boundary file." << endl;

            // Read PtrList of dictionary as dictionary.
            const word oldTypeName = IOPtrList<entry>::typeName;
            const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
            IOPtrList<entry> dictList
            (
                IOobject
                (
                    fieldName,
                    runTime.findInstance
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        fieldName
                    ),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
            // Fake type back to what was in field
            const_cast<word&>(dictList.type()) = dictList.headerClassName();

            // Temporary convert to dictionary
            dictionary fieldDict;
            forAll(dictList, i)
            {
                fieldDict.add(dictList[i].keyword(), dictList[i].dict());
            }

            Info<< "Loaded dictionary " << fieldName
                << " with entries " << fieldDict.toc() << endl;

            // Get the replacement dictionary for the field
            const dictionary& replaceDict = fieldIter().dict();
            Info<< "Merging entries from " << replaceDict.toc() << endl;

            // Merge the replacements in
            merge(fieldDict, replaceDict, literalRE);

            Info<< "fieldDict:" << fieldDict << endl;

            // Convert back into dictList
            wordList doneKeys(dictList.size());

            label nEntries = fieldDict.size();
            forAll(dictList, i)
            {
                doneKeys[i] = dictList[i].keyword();
                dictList.set
                (
                    i,
                    fieldDict.lookupEntry
                    (
                        doneKeys[i],
                        false,
                        true
                    ).clone()
                );
                fieldDict.remove(doneKeys[i]);
            }
            // Add remaining entries
            label sz = dictList.size();
            dictList.setSize(nEntries);
            forAllConstIter(dictionary, fieldDict, iter)
            {
                dictList.set(sz, iter().clone());
            }

            Info<< "Writing modified fieldDict " << fieldName << endl;
            dictList.write();
        }
        else
        {
            // Read dictionary. (disable class type checking so we can load
            // field)
            Info<< "Loading dictionary " << fieldName << endl;
            const word oldTypeName = IOdictionary::typeName;
            const_cast<word&>(IOdictionary::typeName) = word::null;

            IOdictionary fieldDict
            (
                IOobject
                (
                    fieldName,
                    instance,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            const_cast<word&>(IOdictionary::typeName) = oldTypeName;
            // Fake type back to what was in field
            const_cast<word&>(fieldDict.type()) = fieldDict.headerClassName();

            Info<< "Loaded dictionary " << fieldName
                << " with entries " << fieldDict.toc() << endl;

            // Get the replacement dictionary for the field
            const dictionary& replaceDict = fieldIter().dict();
            Info<< "Merging entries from " << replaceDict.toc() << endl;

            // Merge the replacements in
            merge(fieldDict, replaceDict, literalRE);

            Info<< "Writing modified fieldDict " << fieldName << endl;
            fieldDict.regIOobject::write();
        }
    }

    Info<< endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
