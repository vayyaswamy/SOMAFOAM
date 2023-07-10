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

	surfaceScalarField Jif = fvc::interpolate(Ni*Ui) & mesh_.Sf();

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
    	fvm::ddt(Ni, Ui)
   		+ pConvection->fvmDiv(Jif, Ui)
		+ fvm::Sp((Ni*collFrequency_[i]), Ui)
		==
		- fvc::grad(plasmaConstants::boltzC*Ti*Ni)/W(i)
    	+ z_[i]*plasmaConstants::eCharge*Ni*E/W(i)
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
	//Info << "entering updateFlux" << endl;
	volVectorField& Fi = F_[i];

	//Info <<"i = " << i << endl;

	Fi = mu_[i]*(z_[i]*E - plasmaConstants::KBE*fvc::grad(T_[i]));

	//Fi = mu_[i]*(z_[i]*E);

	//Fi.correctBoundaryConditions();

	//Info << "inside " << endl;

	//Info << "exiting updateFlux " << endl;
}


template<class ThermoType>
inline void Foam::mixed<ThermoType>::transportCoeffInterpolate
(
    const label i
)
{


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

	if(transportModel_[i] != "momentum")
	{
		if(diffusionModel_[i] == "eTemp")
		{		
			//Info << "About to do D " << endl;
			D_[i].field() = (1/NG_)*plasmaInterpolateXY
			(
				thermo_.Te().field(),
				graphDiffData_[i].x(),
				graphDiffData_[i].y()
			);
			//Info << "first" << endl;
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

		//Info << "Diffusion interpolated " << endl;
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
									IOobject::NO_READ,
									IOobject::NO_WRITE
								),
								mesh_,
								dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
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
								IOobject::NO_READ,
								IOobject::NO_WRITE
							),
							mesh_,
							dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
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
	

	//Info << "Step 1" << endl;

	if(eonCalculation)
	{	
		EON = mag(E)*1E21/NG_;

		EON.correctBoundaryConditions();
	}

	scalar initialResidual = 1.0;
	
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

				volScalarField& Ni = N_[i] ;

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
						tmp<fvScalarMatrix> NEqn
						(   
							fvm::ddt(Ni)
					   		+ fvm::div((fvc::interpolate(F_[i]) & mesh_.Sf()), Ni, "div(F,Ni)")
							- fvm::laplacian(D_[i], Ni, "laplacian(D,Ni)")
							- fvc::laplacian((mu_[i]*plasmaConstants::KBE*Ni), Ti, "laplacian(D,T)")
							+ fvm::SuSp((-Sy_[i]*plasmaConstants::A/W(i)/Ni), Ni)
						);	
						NEqn->relax();
						NEqn->solve(mesh_.solutionDict().solver("Ni"));
						N_[i].max(1e6);

						yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

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

						tmp<fvScalarMatrix> NEqn
						(   
							fvm::ddt(Ni)
							+ fvm::div((fvc::interpolate(U_[i]) & mesh_.Sf()),Ni, "div(U,Ni)")
							+ fvm::SuSp((-Sy_[i]*plasmaConstants::A/W(i)/Ni), Ni)
						);

						NEqn->relax();

						NEqn->solve(mesh_.solutionDict().solver("Ni"));

						N_[i].max(1e6);

						yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

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

						tmp<fvScalarMatrix> NnEqn
						(   
							fvm::ddt(Ni)
							- fvm::laplacian(D_[i],Ni, "laplacian(D,Nin)")
							+ fvm::SuSp((-Sy_[i]*plasmaConstants::A/W(i)/Ni), Ni)
						);

						NnEqn->relax();
						NnEqn->solve(mesh_.solutionDict().solver("Nin"));
						N_[i].max(1e6);

						yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;
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

		lduSolverPerformance solverPerf;
	
		

		forAll(species(), i)
		{  
			
			//scalar relativeResidual = 1.0;
			if (i != bIndex_ && speciesSolution_[i])
			{
				volScalarField& yi = thermo_.composition().Y(i);

				volScalarField& Ni = N_[i] ;

				

				if(transportModel_[i] != "zeroD")
				{
					transportCoeffInterpolate(i);
				}

				if (i < activeSpecies_)
				{
					if( transportModel_[i] == "driftDiffusion")
					{
						updateFlux(i, E);
						const volScalarField& Ti = T_[i];
						int icorr = 0;

						N_[i].storePrevIter();

						initialResidual = 1.0;
						
						while ((initialResidual >= 1e-5) && (icorr++ <= 3))
						{
						
							tmp<fvScalarMatrix> NEqn
							(   
								fvm::ddt(Ni)
					   			+ fvm::div((fvc::interpolate(F_[i]) & mesh_.Sf()), Ni, "div(F,Ni)")
								- fvm::laplacian(D_[i], Ni, "laplacian(D,Ni)")
								- chemistry.RR(i)*plasmaConstants::A/W(i)
								+ chemistry.dRRDi(i)*Ni
								- fvm::Sp(chemistry.dRRDi(i),Ni)
							);

							//Info << chemistry.dRRDi(i) << endl;
 
							solverPerf = NEqn->solve(mesh_.solutionDict().solver("Ni"));

							initialResidual = solverPerf.initialResidual();

							//Info << "residual = " << initialResidual << endl;

							//N_[i].max(1e4); // setting min value for active species number density

							yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

							updateChemistryCollFreq(chemistry);

						}

						N_[i].relax(); // field relax

						N_[i].max(1e4);

						yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

						updateChemistryCollFreq(chemistry);

						if(diffusionModel_[i] == "einsteinRelation")
						{
							J_[i] == (-mu_[i]*(Ti*fvc::grad(plasmaConstants::KBE*N_[i]) + plasmaConstants::KBE*N_[i]*fvc::grad(Ti) - N_[i]*sign(z_[i])*E));
						}
						else
						{
							J_[i] == (-D_[i]*fvc::grad(N_[i]) - plasmaConstants::KBE*mu_[i]*N_[i]*fvc::grad(Ti) + N_[i]*mu_[i]*sign(z_[i])*E);
						}
						
					}
					else if( transportModel_[i] == "momentum")
					{
						updateVelocity(chemistry, E, i);


						tmp<fvScalarMatrix> NEqn
						(   
							fvm::ddt(Ni)
							+ fvm::div((fvc::interpolate(U_[i]) & mesh_.Sf()),Ni, "div(U,Ni)")
							+ fvm::SuSp((-Sy_[i]*plasmaConstants::A/W(i)/Ni), Ni)
						);

						NEqn->solve(mesh_.solutionDict().solver("Ni"));

						N_[i].max(1e6);

						yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

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

						N_[i].max(1e6);
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


						int icorr = 0;

						initialResidual = 1.0;

						N_[i].storePrevIter();
						
						while ((initialResidual >= 1e-5) && (icorr++ <= 1))
						{
							tmp<fvScalarMatrix> NnEqn
							(   
								fvm::ddt(Ni)
							  - fvm::laplacian(D_[i],Ni, "laplacian(D,Nin)")
							  - chemistry.RR(i)*plasmaConstants::A/W(i)
							  + chemistry.dRRDi(i)*Ni
							  - fvm::Sp(chemistry.dRRDi(i),Ni)
							);

							//Info << "icorr = " << icorr << endl;

							//Info << Ni << endl;

							solverPerf = NnEqn->solve(mesh_.solutionDict().solver("Nin"));

							//Info << "here 1" << endl;

							initialResidual = solverPerf.initialResidual();

							//Info << "here 2" << endl;

							N_[i].max(1e4);

							//Info << "here 3" << endl;

							yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

							//Info << "here 4" << endl;

							updateChemistryCollFreq(chemistry);

							//Info << "here 5" << endl;

						
						}

						N_[i].relax();

						//Info << "here 6" << endl;

						N_[i].max(1e4);

						yi = N_[i]*W(i)/thermo_.rho()/plasmaConstants::A;

						updateChemistryCollFreq(chemistry);
					}
				}
		        yt += yi;  
			}
		}


		volScalarField& yBgas = thermo_.composition().Y()[bIndex_];

		yBgas = scalar(1.0) - yt;

		N_[bIndex_] == thermo_.rho()*thermo_.composition().Y(bIndex_)*plasmaConstants::A/W(bIndex_);

		updateChemistryCollFreq(chemistry);
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
