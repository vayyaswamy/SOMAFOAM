/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "driftDiffusion.H"
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
inline void Foam::driftDiffusion<ThermoType>::updateFlux
(
	const label i,
	const volVectorField& E
)
{     
	volVectorField& Fi = F_[i];

	Fi = thermo_.rho()*mu_[i]*z_[i]*E;

	Fi.correctBoundaryConditions();

	if (restartcapabale && runTime_.write())
	{
		Fi.write();
	}
}

template<class ThermoType>
inline void Foam::driftDiffusion<ThermoType>::transportCoeffInterpolate
(
    const label i
)
{
	if(diffusionModel_[i] == "eTemp")
	{		
		D_[i].field() = (1/NG_)*plasmaInterpolateXY
		(
			thermo_.Te().field(),
			graphDiffData_[i].x(),
			graphDiffData_[i].y()
		);

		D_[i].correctBoundaryConditions();
	}
	else if(diffusionModel_[i] == "EON")
	{		
		D_[i].field() = (1/NG_)*plasmaInterpolateXY
		(
			EON.field(),
			graphDiffData_[i].x(),
			graphDiffData_[i].y()
		);

		D_[i].correctBoundaryConditions();
	}

    if ( i < activeSpecies_)
	{
		if (mobilityModel_[i] == "eTemp")
		{	
			mu_[i].field() = (1/NG_)*plasmaInterpolateXY
			(
				thermo_.Te().field(),
				graphMuData_[i].x(),
				graphMuData_[i].y()
			);

			forAll(mu_[i].boundaryField(), patchI)
			{
				mu_[i].boundaryField()[patchI] = (1/NG_.boundaryField()[patchI])*plasmaInterpolateXY
				(
					thermo_.Te().boundaryField()[patchI],
					graphMuData_[i].x(),
					graphMuData_[i].y()
				);
			}
		}
		else if(mobilityModel_[i] == "EON")
		{	
			mu_[i].field() = (1/NG_)*plasmaInterpolateXY
			(
				EON.field(),
				graphMuData_[i].x(),
				graphMuData_[i].y()
			);

			forAll(mu_[i].boundaryField(), patchI)
			{
				mu_[i].boundaryField()[patchI] = (1/NG_.boundaryField()[patchI])*plasmaInterpolateXY
				(
					EON.boundaryField()[patchI],
					graphMuData_[i].x(),
					graphMuData_[i].y()
				);
			}
		}
		if(diffusionModel_[i] == "einsteinRelation")
		{		
			D_[i] = mu_[i]*plasmaConstants::KBE*T_[i];

			D_[i].correctBoundaryConditions();
		}
	}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::driftDiffusion<ThermoType>::driftDiffusion
(
    hsCombustionThermo& thermo
)
:
    multiSpeciesPlasmaModel(thermo),

    speciesThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),

    EON
    (
		IOobject
		(
			"EON",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    )
{    
    F_.setSize(activeSpecies_);

	updateTemperature();
    
    forAll(F_, i)
    {
        IOobject header
        (
			"F_" + species()[i],
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        );

        if (header.headerOk())
        {
		    F_.set
		    (
		        i, new volVectorField
		        (
		            IOobject
		            (
		                "F_" + species()[i],
		                mesh_.time().timeName(),
		                mesh_,
		                IOobject::MUST_READ,
		                IOobject::NO_WRITE
		            ),
		            mesh_
		        )
		    );
		}
		else
		{
            volVectorField Fdefault
            (
                IOobject
                (
                    "Fdefault",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

		    F_.set
		    (
		        i, new volVectorField
		        (
		            IOobject
		            (
		                "F_" + species()[i],
		                mesh_.time().timeName(),
		                mesh_,
		                IOobject::NO_READ,
		                IOobject::NO_WRITE
		            ),
		            Fdefault
		        )
		    );
		}
    } 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::scalar Foam::driftDiffusion<ThermoType>::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{
	updateTemperature();

	updateChemistryCollFreq(chemistry);

	if(eonCalculation)
	{	
		EON = mag(E)*1E21/NG_;

		EON.correctBoundaryConditions();
	}
	
    volScalarField yt = 0.0*thermo_.composition().Y(0);

	if (multiTimeStep)
	{
		forAll(species(), i)
		{  
			if (i != bIndex_ && speciesSolution_[i])
			{
				if ((timeS_[i] == "LTS") && (fmod(runTime_.deltaT().value(),LTScounter) < SMALL))
				{
					LTSset();
				}

				else if ((timeS_[i] == "MTS") && (fmod(runTime_.deltaT().value(),MTScounter) < SMALL))
				{
					MTSset();
				}

				volScalarField& yi = thermo_.composition().Y(i);

				transportCoeffInterpolate(i);

				if (i < activeSpecies_)
				{
					updateFlux(i, E);

					const volScalarField& Ti = T_[i];

					tmp<fvScalarMatrix> yEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				   		+ fvm::div((fvc::interpolate(F_[i]) & mesh_.Sf()), yi, "div(F,Yi)")
						- fvm::laplacian((thermo_.rho()*D_[i]), yi, "laplacian(D,Yi)")
						- fvc::laplacian((thermo_.rho()*mu_[i]*plasmaConstants::KBE*yi), Ti, "laplacian(D,T)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					yEqn->relax();

					yEqn->solve(mesh_.solutionDict().solver("Yi"));


					yi.max(1e-25);

					N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

					if(diffusionModel_[i] == "einsteinRelation")
					{
						J_[i] == (-mu_[i]*(Ti*fvc::grad(plasmaConstants::KBE*N_[i]) + plasmaConstants::KBE*N_[i]*fvc::grad(Ti) - N_[i]*z_[i]*E));
					}
					else
					{
						J_[i] == (-D_[i]*fvc::grad(N_[i]) - plasmaConstants::KBE*mu_[i]*N_[i]*fvc::grad(Ti) + N_[i]*mu_[i]*z_[i]*E);
					}
				}
				else
				{
					tmp<fvScalarMatrix> ynEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				  		- fvm::laplacian((D_[i]*thermo_.rho()),yi, "laplacian(D,Yin)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					ynEqn->relax();

					ynEqn->solve(mesh_.solutionDict().solver("Yin"));

					yi.max(1e-25);

					N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);
				}
		        yt += yi;  
			}
			STSset();
		}
		volScalarField& yBgas = thermo_.composition().Y()[bIndex_];

		yBgas = scalar(1.0) - yt;

		N_[bIndex_] == thermo_.rho()*thermo_.composition().Y(bIndex_)*plasmaConstants::A/W(bIndex_);
	}
	else
	{
		forAll(species(), i)
		{  
			if (i != bIndex_ && speciesSolution_[i])
			{
				volScalarField& yi = thermo_.composition().Y(i);

				transportCoeffInterpolate(i);

				if (i < activeSpecies_)
				{
					updateFlux(i, E);

					const volScalarField& Ti = T_[i];

					tmp<fvScalarMatrix> yEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				   		+ fvm::div((fvc::interpolate(F_[i]) & mesh_.Sf()), yi, "div(F,Yi)")
						- fvm::laplacian((thermo_.rho()*D_[i]), yi, "laplacian(D,Yi)")
						- fvc::laplacian((thermo_.rho()*mu_[i]*plasmaConstants::KBE*yi), Ti, "laplacian(D,T)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					yEqn->relax();

					yEqn->solve(mesh_.solutionDict().solver("Yi"));

					yi.max(1e-25);

					N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

					if(diffusionModel_[i] == "einsteinRelation")
					{
						J_[i] == (-mu_[i]*(Ti*fvc::grad(plasmaConstants::KBE*N_[i]) + plasmaConstants::KBE*N_[i]*fvc::grad(Ti) - N_[i]*z_[i]*E));
					}
					else
					{
						J_[i] == (-D_[i]*fvc::grad(N_[i]) - plasmaConstants::KBE*mu_[i]*N_[i]*fvc::grad(Ti) + N_[i]*mu_[i]*z_[i]*E);
					}
				}
				else
				{
					tmp<fvScalarMatrix> ynEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				  		- fvm::laplacian((D_[i]*thermo_.rho()),yi, "laplacian(D,Yin)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					ynEqn->relax();

					ynEqn->solve(mesh_.solutionDict().solver("Yin"));

					yi.max(1e-25);

					N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);
				}
		        yt += yi;  
			}
		}

		volScalarField& yBgas = thermo_.composition().Y()[bIndex_];

		yBgas = scalar(1.0) - yt;

		N_[bIndex_] == thermo_.rho()*thermo_.composition().Y(bIndex_)*plasmaConstants::A/W(bIndex_);
	}
    return 0;
}

template<class ThermoType>
bool Foam::driftDiffusion<ThermoType>::read()
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
