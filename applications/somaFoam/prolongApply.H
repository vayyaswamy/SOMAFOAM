// Adaptive mesh refinement, part 2 (after the mesh update): transfer of the
// densities and the electron temperature to refined cells with the
// gradient of the cell they came from (see prolongStore.H).
//
// Densities: n = n0 exp(grad(ln n) . d), with d the displacement from the
// centre of the original cell, scaled so that the cells of one original
// cell keep its content. Electron temperature: Te = Te0 + grad(Te) . d.
// Cells that were not refined have d = 0 and are not changed.

if (prolongCentres.valid())
{
    const vectorField d
    (
        mesh.C().internalField() - prolongCentres().internalField()
    );

    const scalarField& V = mesh.V();

    // Cells created by refinement, grouped by the centre of their origin
    const scalar tol = 1e-6*Foam::cbrt(gMin(V));

    labelList family(mesh.nCells(), -1);
    label nFamilies = 0;
    {
        HashTable<label, word> originToFamily;

        forAll(d, cellI)
        {
            if (mag(d[cellI]) > tol)
            {
                const point& c = prolongCentres().internalField()[cellI];

                OStringStream os;
                os  << label(c.x()/tol) << '_' << label(c.y()/tol) << '_'
                    << label(c.z()/tol);

                const word key(os.str());

                HashTable<label, word>::const_iterator iter =
                    originToFamily.find(key);

                if (iter == originToFamily.end())
                {
                    originToFamily.insert(key, nFamilies);
                    family[cellI] = nFamilies++;
                }
                else
                {
                    family[cellI] = iter();
                }
            }
        }
    }

    if (nFamilies > 0)
    {
        forAll(prolongGradLogN, i)
        {
            // Limited so that a cell changes by a factor of 3 at most
            scalarField factor
            (
                exp
                (
                    max
                    (
                        min
                        (
                            prolongGradLogN[i].internalField() & d,
                            scalar(1.1)
                        ),
                        scalar(-1.1)
                    )
                )
            );

            // Keep the content of each original cell
            scalarField sumV(nFamilies, 0.0);
            scalarField sumVf(nFamilies, 0.0);

            forAll(family, cellI)
            {
                if (family[cellI] >= 0)
                {
                    sumV[family[cellI]] += V[cellI];
                    sumVf[family[cellI]] += V[cellI]*factor[cellI];
                }
            }

            forAll(family, cellI)
            {
                if (family[cellI] >= 0)
                {
                    factor[cellI] *=
                        sumV[family[cellI]]/sumVf[family[cellI]];
                }
                else
                {
                    factor[cellI] = 1;
                }
            }

            volScalarField& Ni = const_cast<volScalarField&>(mspm().N(i));

            Ni.internalField() *= factor;
            Ni.oldTime().internalField() *= factor;
            Ni.correctBoundaryConditions();

            Y[i].internalField() *= factor;
            Y[i].oldTime().internalField() *= factor;
            Y[i].correctBoundaryConditions();
        }

        volScalarField& TeP = thermo.Te();

        scalarField dTe(prolongGradTe().internalField() & d);

        forAll(dTe, cellI)
        {
            // At most 30 % of the value
            const scalar limit = 0.3*TeP[cellI];
            dTe[cellI] = max(min(dTe[cellI], limit), -limit);
        }

        TeP.internalField() += dTe;
        TeP.oldTime().internalField() += dTe;
        TeP.correctBoundaryConditions();
    }
}
