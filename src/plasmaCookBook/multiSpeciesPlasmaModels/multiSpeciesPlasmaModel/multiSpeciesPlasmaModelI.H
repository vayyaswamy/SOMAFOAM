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
  
// ************************************************************************* //
