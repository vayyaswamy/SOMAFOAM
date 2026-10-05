/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Application
    setExpressionFields

Description
    Sets the internal values of scalar fields from formulas of the cell
    centre coordinates, for spatially varying initial conditions. The
    formulas are given in system/setExpressionFieldsDict:

    \verbatim
    variables               // optional: constants and helper formulas
    {
        n0      1e14;
        nPeak   2e16;
        x0      1e-3;
        y0      4e-4;
        w       1.5e-4;
        r2      "sqr(x - x0) + sqr(y - y0)";
    }

    fields
    {
        electron    "n0 + nPeak*exp(-r2/(2*sqr(w)))";
        Arp1        "n0 + nPeak*exp(-r2/(2*sqr(w)))";
        Te          "10000 + 5000*x/2e-3";
    }
    \endverbatim

    x, y and z are the coordinates of the cell centres in metres. The usual
    operators and functions are available (+ - * / ^, sqr, sqrt, exp, log,
    sin, cos, tanh, erf, min, max, mag, pos, neg, ...; pi_ and e_), and a
    formula may use the entries of variables, which may themselves be
    formulas. Each field is read from the time directory, its internal
    values are replaced and it is written back; the boundary conditions and
    their values are kept as they are.

    Options: -time / -latestTime (default: the start time) and -region for
    a region other than the default, e.g. a dielectric.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "equationReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    instantList times = timeSelector::select0(runTime, args);

    if (times.size())
    {
        runTime.setTime(times[times.size() - 1], times.size() - 1);
    }

#   include "createNamedMesh.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

    IOdictionary dict
    (
        IOobject
        (
            "setExpressionFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& fieldsDict = dict.subDict("fields");

    // Cell centre coordinates as the variables of the formulas
    const scalarField x(mesh.C().internalField().component(vector::X));
    const scalarField y(mesh.C().internalField().component(vector::Y));
    const scalarField z(mesh.C().internalField().component(vector::Z));

    equationReader eqns(false);

    eqns.scalarSources().addSource(x, "x");
    eqns.scalarSources().addSource(y, "y");
    eqns.scalarSources().addSource(z, "z");

    if (dict.found("variables"))
    {
        eqns.addSource(dict.subDict("variables"));
    }

    forAllConstIter(dictionary, fieldsDict, iter)
    {
        const word fieldName(iter().keyword());

        volScalarField field
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        eqns.readEquation(fieldsDict, fieldName);

        scalarField values(mesh.nCells(), 0.0);

        eqns.evaluateScalarField(values, fieldName);

        field.internalField() = values;

        Info<< fieldName << ": minimum " << gMin(values)
            << ", maximum " << gMax(values)
            << ", average " << gAverage(values) << endl;

        field.write();
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
