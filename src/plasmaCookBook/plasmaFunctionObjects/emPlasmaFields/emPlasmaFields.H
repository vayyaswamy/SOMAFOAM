/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::emPlasmaFields

Description

SourceFiles
    emPlasmaFields.C
    IOemPlasmaFields.H

\*---------------------------------------------------------------------------*/

#ifndef emPlasmaFields_H
#define emPlasmaFields_H

#include "pointFieldFwd.H"
#include "plasmaConstants.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class emPlasmaFields Declaration
\*---------------------------------------------------------------------------*/

class emPlasmaFields
{
    // Private data

        //- Name of this set of emPlasmaFields objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        wordList objectNames_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        emPlasmaFields(const emPlasmaFields&);

        //- Disallow default bitwise assignment
        void operator=(const emPlasmaFields&);


public:

    //- Runtime type information
    TypeName("emPlasmaFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        emPlasmaFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~emPlasmaFields();


    // Member Functions

        //- Return name of the set of emPlasmaFields
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the emPlasmaFields data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the emPlasmaFields and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
