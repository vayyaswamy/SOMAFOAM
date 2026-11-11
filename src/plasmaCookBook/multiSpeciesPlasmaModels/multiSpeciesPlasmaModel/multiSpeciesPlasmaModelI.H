/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

const inline Foam::speciesTable&
Foam::multiSpeciesPlasmaModel::species() const
{
    return thermo_.composition().species();
}

inline void Foam::multiSpeciesPlasmaModel::LTSset()
{
    runTime_.setDeltaT
    (
        LTScounter
    );

    Info<< "deltaT = " <<  runTime_.deltaT().value() << endl;
}

inline void Foam::multiSpeciesPlasmaModel::MTSset()
{
    runTime_.setDeltaT
    (
        MTScounter
    );
}

inline void Foam::multiSpeciesPlasmaModel::STSset()
{
    runTime_.setDeltaT
    (
        STScounter
    );
}

inline const Foam::volScalarField&
Foam::multiSpeciesPlasmaModel::N(const label i) const
{
    return N_[i];
}

inline const Foam::volScalarField&
Foam::multiSpeciesPlasmaModel::N(const word& specieName) const
{
    return N_[species()[specieName]];
}

inline const Foam::volVectorField&
Foam::multiSpeciesPlasmaModel::U(const label i) const
{
    return U_[i];
}

inline const Foam::volVectorField&
Foam::multiSpeciesPlasmaModel::J(const label i) const
{
    return J_[i];
}

inline const Foam::volVectorField&
Foam::multiSpeciesPlasmaModel::J(const word& specieName) const
{
    return J_[species()[specieName]];
}

inline const Foam::volScalarField&
Foam::multiSpeciesPlasmaModel::mu(const label i) const
{
    return mu_[i];
}

inline void Foam::multiSpeciesPlasmaModel::updateTemperature()
{
    forAll(thermo_.composition().Y(), i)
    {
        const scalarField& TeCells = thermo_.Te().internalField();
        const scalarField& TionCells = thermo_.Tion().internalField();

        scalarField& TiCells = T_[i];

        forAll(TiCells, cellI)
        {
            if (i == eIndex_)
                {
                    TiCells[cellI] = TeCells[cellI];
                }
            else if (i != eIndex_ && i < activeSpecies_)
                {
                    TiCells[cellI] = TionCells[cellI];
                }
        }

        forAll(T_[i].boundaryField(), patchI)
        {
            const fvPatchScalarField& pTe =
                thermo_.Te().boundaryField()[patchI];

            const fvPatchScalarField& pTi =
                thermo_.Tion().boundaryField()[patchI];

            fvPatchScalarField& pTT = T_[i].boundaryField()[patchI];

            forAll(pTT, faceI)
            {
                if (i == eIndex_)
                {
                    pTT[faceI] = pTe[faceI];
                }

                else if (i != eIndex_ && i < activeSpecies_)
                {
                    pTT[faceI] = pTi[faceI];
                }
            }
        }
    }
}

inline void Foam::multiSpeciesPlasmaModel::updateChemistryCollFreq
(
    psiChemistryModel& chemistry
)
{
    if (plasmaChemistry == "temporal")
    {
        // The ODE solve gives only the species rates averaged over the step.
        // The electron energy source and the collision frequencies of the
        // present state come from calculateElectronEnergySource; it is
        // called first because on its first call it falls back to
        // calculateWcf, whose instantaneous rates the solve then replaces
        chemistry.calculateElectronEnergySource();

        chemistry.solve
        (
            runTime_.value() - runTime_.deltaT().value(),
            runTime_.deltaT().value()
        );
    }
    else if (plasmaChemistry == "directSolution")
    {
        chemistry.calculateWcf();
    }
}

inline void Foam::multiSpeciesPlasmaModel::updateElectronEnergyChemistry
(
    psiChemistryModel& chemistry
)
{
    // The passes of the electron energy equation change Te only: the
    // species rates are not needed again
    chemistry.calculateElectronEnergySource();
}

// ************************************************************************* //
