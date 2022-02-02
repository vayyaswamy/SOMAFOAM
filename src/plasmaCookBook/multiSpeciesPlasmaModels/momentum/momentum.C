/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "momentum.H"
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
inline void Foam::momentum<ThermoType>::updateVelocity
(
	const psiChemistryModel& chemistry,
	const volVectorField& E,
	const label i
)
{     
	const volScalarField& Ti = T_[i];

	const volScalarField& Ni = N_[i];

	volScalarField Yi = thermo_.rho()*thermo_.composition().Y(i);

	volVectorField& Ui = U_[i];

	surfaceScalarField Jif = fvc::interpolate(Yi*Ui) & mesh_.Sf();

	if (collisionFrequency_[i] == "rrBased")
	{
		collFrequency_[i] = chemistry.collFreq(i);
    }
	else if (collisionFrequency_[i] == "muBased")
	{
		collFrequency_[i] = plasmaConstants::eChargeA/mu_[i]/W(i);
	}

	tmp<fv::convectionScheme<vector> > pConvection
	(
		fv::convectionScheme<vector>::New
		(
		    mesh_,
			Jif,
		    mesh_.schemesDict().divScheme("div(flux,Ui_h)")
		)
	);

	tmp<fvVectorMatrix> uEqn
	(   
    	fvm::ddt(Yi, Ui)
   		+ pConvection->fvmDiv(Jif, Ui)
		+ fvm::Sp((Yi*collFrequency_[i]), Ui)
		==
		- fvc::grad(plasmaConstants::boltzC*Ti*Ni)
    	+ z_[i]*plasmaConstants::eCharge*Ni*E
	);

	uEqn->relax();

    uEqn->solve(mesh_.solutionDict().solver("Ui"));

	if (restartcapabale && runTime_.write())
	{
		Ui.write();
	}
}

template<class ThermoType>
inline void Foam::momentum<ThermoType>::transportCoeffInterpolate
(
    const label i
)
{
    if (i < activeSpecies_ && collisionFrequency_[i] == "muBased")
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
	}
	else
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
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::momentum<ThermoType>::momentum
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
    U_.setSize(activeSpecies_);

	updateTemperature();
    
    forAll(U_, i)
    {
        IOobject header
        (
			"U_" + species()[i],
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        );

        if (header.headerOk())
        {
		    U_.set
		    (
		        i, new volVectorField
		        (
		            IOobject
		            (
		                "U_" + species()[i],
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
            volVectorField Udefault
            (
                IOobject
                (
                    "Udefault",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

		    U_.set
		    (
		        i, new volVectorField
		        (
		            IOobject
		            (
		                "U_" + species()[i],
		                mesh_.time().timeName(),
		                mesh_,
		                IOobject::NO_READ,
		                IOobject::NO_WRITE
		            ),
		            Udefault
		        )
		    );
		}
    } 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::scalar Foam::momentum<ThermoType>::correct
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
					updateVelocity(chemistry, E, i);

					tmp<fvScalarMatrix> yEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				   		+ fvm::div(fvc::interpolate(U_[i]) & mesh_.Sf(), yi, "div(F,Yi)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					yEqn->relax();

					yEqn->solve(mesh_.solutionDict().solver("Yi"));

					yi.max(0.0);

					N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

					J_[i] = N_[i]*U_[i];
				}
				else
				{
					tmp<fvScalarMatrix> ynEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				  		- fvm::laplacian((D_[i]*thermo_.rho()),yi, "laplacian(D,Yin)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					ynEqn->solve(mesh_.solutionDict().solver("Yin"));

					yi.max(0.0);

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
					updateVelocity(chemistry, E, i);

					tmp<fvScalarMatrix> yEqn
					(   
						fvm::ddt(thermo_.rho(),yi)
				   		+ fvm::div(fvc::interpolate(U_[i]) & mesh_.Sf(), yi, "div(F,Yi)")
						+ fvm::SuSp((-Sy_[i]/yi), yi)
					);

					yEqn->relax();

					yEqn->solve(mesh_.solutionDict().solver("Yi"));

					yi.max(0.0);

					N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

					J_[i] == N_[i]*U_[i];
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

					yi.max(0.0);

					N_[i] = thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);
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
bool Foam::momentum<ThermoType>::read()
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

