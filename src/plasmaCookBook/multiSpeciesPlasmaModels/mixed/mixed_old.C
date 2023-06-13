/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "mixed.H"
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
template<class ThermoType>
void Foam::mixed<ThermoType>::mixedModelInput
()
{
	IOdictionary diffDict
    (
        IOobject
        (
            "plasmaProperties", 
            thermo_.T().mesh().time().constant(),
            thermo_.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE 
        )
    );

    forAll(species(), i)
    {
		const word& currentSpecie = species()[i];
	    const dictionary& subDict = diffDict.subDict(currentSpecie);

        transportModel_.set
        (
             i,
             new word("none")
        );

        transportModel_[i] =
        (
			subDict.lookupOrDefault<word>("transportModel", "none")
        );

		if (i < activeSpecies_)
		{
			if (transportModel_[i] == "momentum")
			{
				momentumCount++;
			}
			else if (transportModel_[i] == "driftDiffusion")
			{
				ddCount++;
			}
		}
	} 
}


template<class ThermoType>
inline void Foam::mixed<ThermoType>::updateVelocity
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
inline void Foam::mixed<ThermoType>::updateFlux
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
inline void Foam::mixed<ThermoType>::transportCoeffInterpolate
(
    const label i
)
{
	if(transportModel_[i] != "momentum")
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
		else if(diffusionModel_[i] == "einsteinRelation")
		{		
			D_[i] == mu_[i]*plasmaConstants::KBE*T_[i];
		}
	}

    if ((i < activeSpecies_ && transportModel_[i] != "momentum") || (transportModel_[i] == "momentum" && collisionFrequency_[i] == "muBased"))
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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
template<class ThermoType>
Foam::mixed<ThermoType>::mixed
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

	momentumCount(0.0),

	ddCount(0.0),

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
    transportModel_.setSize(species().size());

	mixedModelInput();

    U_.setSize(momentumCount);

    F_.setSize(ddCount);

	updateTemperature();
    
	forAll(species(), i)
	{  
		if (i < activeSpecies_)
		{
			if( transportModel_[i] != "zeroD" )
			{
				IOobject headerU
				(
					"U_" + species()[i],
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ
				);

				IOobject headerF
				(
					"F_" + species()[i],
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ
				);

				if( transportModel_[i] == "momentum" )
				{
					if (headerU.headerOk())
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
									IOobject::AUTO_WRITE
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
									IOobject::AUTO_WRITE
								),
								Udefault
							)
						);
					}
				}
				else if( transportModel_[i] == "driftDiffusion" )
				{
					if (headerF.headerOk())
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
		}	
    } 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::scalar Foam::mixed<ThermoType>::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{
	updateTemperature();

	updateChemistryCollFreq(chemistry);

	//Info << "Step 1" << endl;

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

				if( transportModel_[i] != "zeroD")
				{
					transportCoeffInterpolate(i);
				}

				if (i < activeSpecies_)
				{
					if( transportModel_[i] == "driftDiffusion")
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

						yi.max(0.0);

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e6);

						if(diffusionModel_[i] == "einsteinRelation")
						{
							J_[i] == (-mu_[i]*(Ti*fvc::grad(plasmaConstants::KBE*N_[i]) + plasmaConstants::KBE*N_[i]*fvc::grad(Ti) - N_[i]*z_[i]*E));
						}
						else
						{
							J_[i] == (-D_[i]*fvc::grad(N_[i]) - plasmaConstants::KBE*mu_[i]*N_[i]*fvc::grad(Ti) + N_[i]*mu_[i]*z_[i]*E);
						}
					}
					else if( transportModel_[i] == "momentum")
					{
						updateVelocity(chemistry, E, i);

						tmp<fvScalarMatrix> yEqn
						(   
							fvm::ddt(thermo_.rho(),yi)
							+ fvm::div((fvc::interpolate(U_[i]*thermo_.rho()) & mesh_.Sf()),yi, "div(F,Yi)")
							+ fvm::SuSp((-Sy_[i]/yi), yi)
						);

						yEqn->relax();

						yEqn->solve(mesh_.solutionDict().solver("Yi"));

						yi.max(0.0);

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e6);

						J_[i] == N_[i]*U_[i];
					}
					else if( transportModel_[i] == "zeroD")
					{
						tmp<fvScalarMatrix> yEqn
						(   
							fvm::ddt(thermo_.rho(),yi) == Sy_[i]
						);

						yEqn->solve(mesh_.solutionDict().solver("Yi"));

						yi.max(0.0);

						

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);
					}
				}
				else
				{
					if( transportModel_[i] == "zeroD")
					{
						tmp<fvScalarMatrix> ynEqn
						(   
							fvm::ddt(thermo_.rho(),yi) == Sy_[i]
						);

						ynEqn->solve(mesh_.solutionDict().solver("Yin"));

						yi.max(0.0);

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);
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

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);
					}
				}
		        yt += yi;  
			}
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



				if(transportModel_[i] != "zeroD")
				{
					//Info << "Step 2 coming" << endl;
					transportCoeffInterpolate(i);
					//Info << "Step 2 done" << endl;
				}

				//Info << "Step 2 " << endl;

				if (i < activeSpecies_)
				{
					if( transportModel_[i] == "driftDiffusion")
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

						//Info << "Step 3 " << endl;

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);

						if(diffusionModel_[i] == "einsteinRelation")
						{
							J_[i] == (-mu_[i]*(Ti*fvc::grad(plasmaConstants::KBE*N_[i]) + plasmaConstants::KBE*N_[i]*fvc::grad(Ti) - N_[i]*z_[i]*E));
						}
						else
						{
							J_[i] == (-D_[i]*fvc::grad(N_[i]) - plasmaConstants::KBE*mu_[i]*N_[i]*fvc::grad(Ti) + N_[i]*mu_[i]*z_[i]*E);
						}
					}
					else if( transportModel_[i] == "momentum")
					{
						updateVelocity(chemistry, E, i);

						tmp<fvScalarMatrix> yEqn
						(   
							fvm::ddt(thermo_.rho(),yi)
							+ fvm::div((fvc::interpolate(U_[i]*thermo_.rho()) & mesh_.Sf()),yi, "div(F,Yi)")
							+ fvm::SuSp((-Sy_[i]/yi), yi)
						);

						yEqn->relax();

						yEqn->solve(mesh_.solutionDict().solver("Yi"));

						

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);

						J_[i] == N_[i]*U_[i];
					}
					else if( transportModel_[i] == "zeroD")
					{
						tmp<fvScalarMatrix> yEqn
						(   
							fvm::ddt(thermo_.rho(),yi) == Sy_[i]
						);

						yEqn->solve(mesh_.solutionDict().solver("Yi"));

						

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);
					}
				}
				else
				{
					if( transportModel_[i] == "zeroD")
					{
						tmp<fvScalarMatrix> ynEqn
						(   
							fvm::ddt(thermo_.rho(),yi) == Sy_[i]
						);

						ynEqn->solve(mesh_.solutionDict().solver("Yin"));

						yi.max(0.0);

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);
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

						

						N_[i] == thermo_.rho()*(thermo_.composition().Y(i)*plasmaConstants::A)/W(i);

						N_[i].max(1e10);
					}
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
bool Foam::mixed<ThermoType>::read()
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
