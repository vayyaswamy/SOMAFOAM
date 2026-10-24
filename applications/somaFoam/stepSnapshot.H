/*---------------------------------------------------------------------------*\
Class
    Foam::stepSnapshot

Description
    Values of all fields registered on a mesh at the start of a time step,
    so that a rejected step can be undone and retried with a smaller step.

    The internal and boundary values of every registered vol and surface
    scalar and vector field are copied. Old-time fields ("<name>_0") are not:
    a retry stays at the same time index, so they keep the values of the
    start of the step. Fields that no longer exist on restore (temporaries)
    are skipped, and so are the mesh geometry fields.

\*---------------------------------------------------------------------------*/

#ifndef stepSnapshot_H
#define stepSnapshot_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, template<class> class PatchField, class GeoMesh>
class geoFieldSnapshot
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;

    const objectRegistry& db_;

    wordList names_;
    List<Field<Type> > internal_;
    List<List<Field<Type> > > boundary_;

    //- Old-time fields, and the mesh geometry and interpolation fields
    //  (registered as vol and surface fields, some of them slices of the
    //  mesh data), are not state of the step
    static bool skip(const word& name)
    {
        static const char* geometry[] =
        {
            "C", "Cf", "Sf", "magSf", "weights", "deltaCoeffs",
            "nonOrthDeltaCoeffs", "nonOrthCorrectionVectors", "magLongDelta",
            "meshPhi"
        };

        for (unsigned i = 0; i < sizeof(geometry)/sizeof(geometry[0]); i++)
        {
            if (name == geometry[i])
            {
                return true;
            }
        }

        return name.size() > 2 && name.substr(name.size() - 2) == "_0";
    }

public:

    geoFieldSnapshot(const objectRegistry& db)
    :
        db_(db)
    {}

    void store()
    {
        const HashTable<const fieldType*> fields
        (
            db_.lookupClass<fieldType>()
        );

        names_.setSize(fields.size());
        label n = 0;

        forAllConstIter(typename HashTable<const fieldType*>, fields, iter)
        {
            if (!skip(iter.key()))
            {
                names_[n++] = iter.key();
            }
        }

        names_.setSize(n);
        internal_.setSize(n);
        boundary_.setSize(n);

        forAll(names_, i)
        {
            const fieldType& f = *fields[names_[i]];

            internal_[i] = f.internalField();

            boundary_[i].setSize(f.boundaryField().size());

            forAll(f.boundaryField(), patchi)
            {
                boundary_[i][patchi] = f.boundaryField()[patchi];
            }
        }
    }

    void restore() const
    {
        forAll(names_, i)
        {
            if (!db_.foundObject<fieldType>(names_[i]))
            {
                continue;
            }

            fieldType& f =
                const_cast<fieldType&>(db_.lookupObject<fieldType>(names_[i]));

            f.internalField() = internal_[i];

            forAll(f.boundaryField(), patchi)
            {
                f.boundaryField()[patchi] == boundary_[i][patchi];
            }
        }
    }
};


class stepSnapshot
{
    geoFieldSnapshot<scalar, fvPatchField, volMesh> volScalar_;
    geoFieldSnapshot<vector, fvPatchField, volMesh> volVector_;
    geoFieldSnapshot<scalar, fvsPatchField, surfaceMesh> surfaceScalar_;
    geoFieldSnapshot<vector, fvsPatchField, surfaceMesh> surfaceVector_;

public:

    stepSnapshot(const fvMesh& mesh)
    :
        volScalar_(mesh),
        volVector_(mesh),
        surfaceScalar_(mesh),
        surfaceVector_(mesh)
    {}

    void store()
    {
        volScalar_.store();
        volVector_.store();
        surfaceScalar_.store();
        surfaceVector_.store();
    }

    void restore() const
    {
        volScalar_.restore();
        volVector_.restore();
        surfaceScalar_.restore();
        surfaceVector_.restore();
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
