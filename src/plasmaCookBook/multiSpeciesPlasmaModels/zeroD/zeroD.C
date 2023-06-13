/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "zeroD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::zeroD<ThermoType>::zeroD
(
    hsCombustionThermo& thermo
)
:
    multiSpeciesPlasmaModel(thermo),

    speciesThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::scalar Foam::zeroD<ThermoType>::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{
    volScalarField yt = 0.0*thermo_.composition().Y(0);

	updateChemistryCollFreq(chemistry);

	if (multiTimeStep)
	{
		FatalError << "0-D model not implemented for multiple time scales" << nl << abort(FatalError);
	}
	else
	{
		forAll(species(), i)
		{  
			if (i != bIndex_ && speciesSolution_[i])
			{
				volScalarField& yi = thermo_.composition().Y(i);

				tmp<fvScalarMatrix> yEqn
				(   
					fvm::ddt(thermo_.rho(),yi) == Sy_[i]
				);

				yEqn->solve(mesh_.solutionDict().solver("Yi"));
	
				yi.max(0.0);

		        yt += yi;  

				N_[i] == thermo_.rho()*thermo_.composition().Y(i)*plasmaConstants::A/W(i);
			}
		}
		volScalarField& yBgas = thermo_.composition().Y()[bIndex_];

		yBgas == 1 - yt;

		N_[bIndex_] == thermo_.rho()*thermo_.composition().Y(bIndex_)*plasmaConstants::A/W(bIndex_);
	}
    return 0;
}

template<class ThermoType>
bool Foam::zeroD<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
