/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "refinementIndicators.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Size of the cells in the directions that the mesh resolves: the volume
//  divided by the extent in the empty directions, to the power one over the
//  number of resolved directions
tmp<scalarField> cellSize(const fvMesh& mesh)
{
    const Vector<label>& directions = mesh.geometricD();

    const vector span(mesh.bounds().span());

    scalar thickness = 1;

    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        if (directions[cmpt] == -1)
        {
            thickness *= span[cmpt];
        }
    }

    return pow(mesh.V().field()/thickness, 1.0/mesh.nGeometricD());
}


//- Add the contribution of one field to the indicator (maximum is kept)
template<class Type>
void addIndicator
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const word& type,
    const scalar weight,
    const scalar scale,
    const scalar floor,
    scalarField& indicator
)
{
    const Field<Type>& f = fld.internalField();

    if (type == "magnitude")
    {
        forAll(f, cellI)
        {
            indicator[cellI] =
                max(indicator[cellI], weight*mag(f[cellI])/scale);
        }
    }
    else if (type == "relativeGradient")
    {
        // Cell size times the gradient relative to the value: the cell
        // size over the local scale length of the field
        const scalarField magGrad(mag(fvc::grad(fld))().internalField());

        const scalarField h(cellSize(mesh));

        forAll(f, cellI)
        {
            indicator[cellI] = max
            (
                indicator[cellI],
                weight*h[cellI]*magGrad[cellI]/(mag(f[cellI]) + floor + VSMALL)
            );
        }
    }
    else if (type == "gradient")
    {
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        forAll(nei, faceI)
        {
            const label a = own[faceI];
            const label b = nei[faceI];

            const scalar value = weight*mag(f[a] - f[b])/scale;

            indicator[a] = max(indicator[a], value);
            indicator[b] = max(indicator[b], value);
        }
    }
    else
    {
        FatalErrorIn("refinementIndicators::calculate(...)")
            << "Unknown indicator type " << type << " for field "
            << fld.name() << nl << "Valid types are relativeGradient,"
            << " gradient, magnitude and debyeLength" << abort(FatalError);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::refinementIndicators::calculate
(
    const fvMesh& mesh,
    const dictionary& dict,
    scalarField& indicator
)
{
    indicator.setSize(mesh.nCells());
    indicator = 0.0;

    const PtrList<dictionary> entries(dict.lookup("indicators"));

    forAll(entries, i)
    {
        const dictionary& entry = entries[i];

        const word type(entry.lookup("type"));
        const word fieldName(entry.lookup("field"));
        const scalar weight = entry.lookupOrDefault<scalar>("weight", 1.0);
        const scalar scale = entry.lookupOrDefault<scalar>("scale", 1.0);
        const scalar floor = entry.lookupOrDefault<scalar>("floor", 0.0);

        if (type == "debyeLength")
        {
            // Cell size over the electron Debye length
            const word TeName(entry.lookupOrDefault<word>("temperature", "Te"));

            if
            (
                !mesh.foundObject<volScalarField>(fieldName)
             || !mesh.foundObject<volScalarField>(TeName)
            )
            {
                FatalErrorIn("refinementIndicators::calculate(...)")
                    << "debyeLength needs the volScalarFields " << fieldName
                    << " (density) and " << TeName << " (temperature);"
                    << " available: " << mesh.names<volScalarField>()
                    << abort(FatalError);
            }

            const scalarField& n =
                mesh.lookupObject<volScalarField>(fieldName).internalField();
            const scalarField& Te =
                mesh.lookupObject<volScalarField>(TeName).internalField();

            const scalarField h(cellSize(mesh));

            const scalar eps0 = 8.8541878128e-12;
            const scalar kB = 1.380649e-23;
            const scalar qe = 1.602176634e-19;

            forAll(n, cellI)
            {
                if (n[cellI] > 0)
                {
                    // Floor on Te avoids a zero Debye length in cold cells
                    const scalar lambdaD = Foam::sqrt
                    (
                        eps0*kB*max(Te[cellI], 300.0)/(qe*qe*n[cellI])
                    );

                    indicator[cellI] =
                        max(indicator[cellI], weight*h[cellI]/lambdaD);
                }
            }

            continue;
        }

        if (mesh.foundObject<volScalarField>(fieldName))
        {
            addIndicator
            (
                mesh,
                mesh.lookupObject<volScalarField>(fieldName),
                type,
                weight,
                scale,
                floor,
                indicator
            );
        }
        else if (mesh.foundObject<volVectorField>(fieldName))
        {
            addIndicator
            (
                mesh,
                mesh.lookupObject<volVectorField>(fieldName),
                type,
                weight,
                scale,
                floor,
                indicator
            );
        }
        else
        {
            FatalErrorIn("refinementIndicators::calculate(...)")
                << "Indicator field " << fieldName << " not found among the"
                << " volScalarFields " << mesh.names<volScalarField>()
                << " and volVectorFields " << mesh.names<volVectorField>()
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
