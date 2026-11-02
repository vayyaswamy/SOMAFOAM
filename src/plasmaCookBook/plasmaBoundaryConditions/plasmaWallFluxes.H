/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Description
    Quantities at a wall shared by the electron wall conditions
    (driftDiffusionElectronDensity, electronTemperature):

    - the ion flux to the wall that drives secondary emission, computed at
      the wall from the ion densities and the field (emissionIonFlux wall):
          sum_i (mu_i n_i E.n - D_i dn_i/dn), positive part,
      with the wall values of the densities and the outward normal;
    - the electron temperature of the thermal speed of the wall flux, at the
      wall face (thermalSpeed wall) or in the wall cell (thermalSpeed cell).

\*---------------------------------------------------------------------------*/

#ifndef plasmaWallFluxes_H
#define plasmaWallFluxes_H

#include "fvPatchFields.H"
#include "volFields.H"

namespace Foam
{

//- Ion flux to the wall from the ion densities N_<ion> and the field E
inline tmp<scalarField> wallIonFlux(const fvPatch& p, const wordList& ions)
{
    tmp<scalarField> tG(new scalarField(p.size(), 0.0));
    scalarField& G = tG();

    const vectorField n(p.nf());
    const scalarField En(p.lookupPatchField<volVectorField, vector>("E") & n);

    forAll(ions, i)
    {
        const fvPatchScalarField& Ni =
            p.lookupPatchField<volScalarField, scalar>("N_" + ions[i]);
        const fvPatchScalarField& mui =
            p.lookupPatchField<volScalarField, scalar>("mu_" + ions[i]);
        const fvPatchScalarField& Di =
            p.lookupPatchField<volScalarField, scalar>("D_" + ions[i]);

        G += mui*Ni*En - Di*Ni.snGrad();
    }

    G = max(G, scalar(0));

    return tG;
}


//- Electron temperature of the thermal speed of the wall flux
inline tmp<scalarField> wallThermalTemperature
(
    const fvPatchScalarField& Te,
    const bool cell
)
{
    return cell ? Te.patchInternalField() : tmp<scalarField>(new scalarField(Te));
}

} // End namespace Foam

#endif

// ************************************************************************* //
