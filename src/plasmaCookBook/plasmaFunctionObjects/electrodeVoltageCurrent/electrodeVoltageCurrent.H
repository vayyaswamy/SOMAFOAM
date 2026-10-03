/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::electrodeVoltageCurrent

Description
    Writes the voltage and the current of selected electrode patches as a
    function of time, one file per patch:

        <case>/<name>/<start time>/<patch>.dat

    with the columns

        t [s]  V [V]  I_conduction [A]  I_displacement [A]  I_total [A]
        j_total [A/m2]

    V is the area average of the potential on the patch. The currents are
    the patch integrals of the conduction current density (charged-species
    fluxes times charge) and of the total current density (conduction plus
    displacement, eps dE/dt), as computed by the solver; they are positive
    when current flows from the electrode into the plasma. j_total is the
    total current divided by the patch area.

    Example (system/controlDict):
    \verbatim
    functions
    {
        electrodes
        {
            type                electrodeVoltageCurrent;
            functionObjectLibs  ("libplasmaFunctionObjects.so");
            patches             (electrode ground);

            // optional
            outputInterval      1;        // time steps between lines
            region              region0;
            potential           Phi;
            conductionCurrent   Jnet;
            totalCurrent        Jtot;
        }
    }
    \endverbatim

SourceFiles
    electrodeVoltageCurrent.C

\*---------------------------------------------------------------------------*/

#ifndef electrodeVoltageCurrent_H
#define electrodeVoltageCurrent_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "PtrList.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class electrodeVoltageCurrent Declaration
\*---------------------------------------------------------------------------*/

class electrodeVoltageCurrent
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

        //- Electrode patches
        wordList patchNames_;

        //- Field names
        word potentialName_;
        word conductionCurrentName_;
        word totalCurrentName_;

        //- Time steps between output lines
        label outputInterval_;

        //- Output files, one per patch (master processor only)
        PtrList<OFstream> files_;


    // Private Member Functions

        //- Open the output files and write the headers
        void makeFiles();

        //- Write one line per patch
        void write();

        //- Disallow default bitwise copy construct
        electrodeVoltageCurrent(const electrodeVoltageCurrent&);

        //- Disallow default bitwise assignment
        void operator=(const electrodeVoltageCurrent&);


public:

    //- Runtime type information
    TypeName("electrodeVoltageCurrent");


    // Constructors

        //- Construct from components
        electrodeVoltageCurrent
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
