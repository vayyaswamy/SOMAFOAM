/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "multiSpeciesPlasmaModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(multiSpeciesPlasmaModel, 0);
    defineRunTimeSelectionTable(multiSpeciesPlasmaModel, fvMesh);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //  
void Foam::multiSpeciesPlasmaModel::updateChemistry
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
		chemistry.calculate();
	}
}
     
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSpeciesPlasmaModel::multiSpeciesPlasmaModel
(
    hsCombustionThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            "plasmaProperties",
			thermo.T().mesh().time().constant(),
			thermo.T().mesh(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
        )
    ),

    thermo_(thermo),

    mesh_(thermo.T().mesh()),

	runTime_(const_cast<Time&>(mesh_.time())),

    activeSpecies_(readInt(IOdictionary::lookup("activeSpecies"))),

	multiTimeStep(IOdictionary::lookupOrDefault("multiTimeStep", false)),

	STScounter(runTime_.deltaT().value()),

	MTScounter(IOdictionary::lookupOrDefault<scalar>("LTSvalue", 0.0)),

	LTScounter(IOdictionary::lookupOrDefault<scalar>("MTSvalue", 0.0)),

	eSpecie("electron"),

	eIndex_(species()[eSpecie]),

	bGas(IOdictionary::lookup("backgroundGas")),

	bIndex_(species()[bGas]),

	plasmaChemistry("directSolution"),

	restartcapabale(runTime_.controlDict().lookup("restartCapable")),

	eonCalculation(false),

    NG_
    (
		IOobject
		(
			"NG_",
            mesh_.time().timeName(),
            mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    )
{
    graphDiffData_.setSize(species().size());

    graphMuData_.setSize(activeSpecies_);

    z_.setSize(activeSpecies_);

    diffusionModel_.setSize(species().size());

    mobilityModel_.setSize(activeSpecies_);

	collisionFrequency_.setSize(activeSpecies_);

	transportModel_.setSize(species().size());

	speciesSolution_.setSize(species().size());

    timeS_.setSize(species().size());

    Sy_.setSize(species().size());

    D_.setSize(species().size());

    mu_.setSize(activeSpecies_);

    T_.setSize(species().size());

    J_.setSize(activeSpecies_);

    N_.setSize(species().size());

    collFrequency_.setSize(activeSpecies_);
      
    forAll(thermo.composition().Y(), i)
    {
        Sy_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "Sy_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0),
	            zeroGradientFvPatchScalarField::typeName
            )
        );

        D_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "D_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0),
	            zeroGradientFvPatchScalarField::typeName
            )
        );

        //Info << "species " << species()[i] << endl;

        N_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    //"N_" + species()[i],
                    species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
                //dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
            )
        );

        N_[i].rename("N_" + species()[i]);

        //Info << "N = " << N_[i] << endl;

        T_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "T_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.0)
            )
        );
		
    	if ( i < activeSpecies_)
		{
		    mu_.set
		    (
		        i, new volScalarField
		        (
		            IOobject
		            (
		                "mu_" + species()[i],
		                mesh_.time().timeName(),
		                mesh_,
		                IOobject::NO_READ,
		                IOobject::NO_WRITE
		            ),
		            mesh_,
		            dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
		        )
		    );

		    collFrequency_.set
		    (
		        i, new volScalarField
		        (
		            IOobject
		            (
		                "collFrequency_" + species()[i],
		                mesh_.time().timeName(),
		                mesh_,
		                IOobject::NO_READ,
		                IOobject::NO_WRITE
		            ),
		            mesh_,
		            dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
		        )
		    );

		    J_.set
		    (
		        i, new volVectorField
		        (
		            IOobject
		            (
		                "J_" + species()[i],
		                mesh_.time().timeName(),
		                mesh_,
		                IOobject::NO_READ,
		                IOobject::NO_WRITE
		            ),
		            mesh_,
		            dimensionedVector("zero", dimensionSet(1, -1, -1, 0, 0), vector::zero)
		        )
		    );
		}        
    } 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::multiSpeciesPlasmaModel::input
()
{ 
	IOdictionary diffDict
    (
        IOobject
        (
            "plasmaProperties", 
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE  
        )
    );

	IOdictionary chemDict
    (
        IOobject
        (
            "chemistryProperties", 
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE  
        )
    );

	plasmaChemistry = (chemDict.lookupOrDefault<word>("plasmaChemistryModel", "directSolution"));

	if ((plasmaChemistry != "temporal") && (plasmaChemistry != "directSolution"))
	{
		WarningIn
		(
		    "plasmaChemistry"
		)   << "in chemistryProperties, unknown plasmaChemistry type."
		    << endl;
	}

	int a(0);

    forAll(species(), i)
    {
		const word& currentSpecie = species()[i];
        const dictionary& subDict = diffDict.subDict(currentSpecie);

        diffusionModel_.set
        (
             i,
             new word("none")
        );

        diffusionModel_[i] =
        (
			subDict.lookupOrDefault<word>("diffusionModel", "none")
        );

		if ((diffusionModel_[i] != "EON") && (diffusionModel_[i] != "eTemp") && (diffusionModel_[i] != "constant") && (diffusionModel_[i] != "einsteinRelation") && (diffusionModel_[i] != "none"))
		{
			WarningIn
			(
			    "diffusionModel"
			)   << "unknown diffusionModel type."
			    << endl;
		}

        speciesSolution_.set
        (
             i,
             new Switch(subDict.lookup("speciesSolution"))
        );

		if (multiTimeStep)
		{
		    timeS_.set
		    (
		         i,
		         new word(subDict.lookup("timeScale"))
		    );
		}

		if (diffusionModel_[i] == "constant")
		{
			tmp<volScalarField> d_
			(
				new volScalarField
				(
					IOobject
					(
					    "d_",
			            mesh_.time().timeName(),
			            mesh_,
			            IOobject::NO_READ,
			            IOobject::NO_WRITE
					),
					mesh_,
			        dimensionedScalar("d_", dimensionSet(1, -1, -1, 0, 0, 0, 0), subDict.lookup("D"))
				)
			);

			D_[i] = d_;

			D_[i].correctBoundaryConditions();
		}
		else if (diffusionModel_[i] == "eTemp" || diffusionModel_[i] == "EON")
		{
			IFstream file_D(thermo_.T().mesh().time().constant()/"D_" + species()[i]);

			graphDiffData_.set(i, new graph("D_data_file","inter_data","D_data", file_D));;
		}

		//N_[i] = thermo_.rho()*thermo_.composition().Y(i)*plasmaConstants::A/W(i);

		//N_[i].correctBoundaryConditions();


        if ( i < activeSpecies_)
		{
		    z_.set
		    (
		         i,
		         new scalar(0.0)
		    );

		    z_[i] = (subDict.lookupOrDefault<scalar>("charge", 0.0));

		    collisionFrequency_.set
		    (
		         i,
		         new word("none")
		    );

			collisionFrequency_[i] = (subDict.lookupOrDefault<word>("collisionFrequency", "none"));

			if ((collisionFrequency_[i] != "muBased") && (collisionFrequency_[i] != "rrBased") && (collisionFrequency_[i] != "none"))
			{
				WarningIn
				(
				    "collisionFrequency"
				)   << "unknown collisionFrequency type."
				    << endl;
			}

		    mobilityModel_.set
		    (
		         i,
		         new word("none")
		    );

		    mobilityModel_[i] =
		    (
				subDict.lookupOrDefault<word>("mobilityModel", "none")
		    );

			if ((mobilityModel_[i] != "EON") && (mobilityModel_[i] != "eTemp") && (mobilityModel_[i] != "constant") && (mobilityModel_[i] != "none"))
			{
				WarningIn
				(
				    "mobilityModel"
				)   << "unknown mobilityModel type."
				    << endl;
			}

			if (mobilityModel_[i] == "constant")
			{
				tmp<volScalarField> mui_
				(
					new volScalarField
					(
						IOobject
						(
							"mui_",
					        mesh_.time().timeName(),
					        mesh_,
					        IOobject::NO_READ,
					        IOobject::NO_WRITE
						),
						mesh_,
					    dimensionedScalar("mui_", dimensionSet(1, -1, -1, 0, 0, 0, 0), subDict.lookup("mu"))
					)
				);

				forAll(mu_[i], celli)
				{
					mu_[i][celli] = mui_()[celli];
				}

				forAll(mu_[i].boundaryField(), patchi)
				{
					const fvPatchScalarField& pmu =
						mui_().boundaryField()[patchi];

					fvPatchScalarField& pmui = mu_[i].boundaryField()[patchi];

					forAll(pmui, facei)
					{
						pmui[facei] = pmu[facei];
					}
				}
			}
			else if (mobilityModel_[i] == "eTemp" || mobilityModel_[i] == "EON")
			{
				IFstream file_mu(thermo_.T().mesh().time().constant()/"mu_" + species()[i]);

				graphMuData_.set(i, new graph("mu_data_file","inter_data","mu_data", file_mu));
			}
		}

		if((diffusionModel_[i].find("EON")) || (mobilityModel_[i].find("EON"))) 
		{
			a++; 
		}
    }
	
	if(a!=0)
	{
		eonCalculation=true;
	}
}

Foam::scalar Foam::multiSpeciesPlasmaModel::correct
(
	PtrList<volScalarField>& Y,
    const volVectorField& E,
    psiChemistryModel& chemistry,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{ 
	updateTemperature();

	updateChemistryCollFreq(chemistry);
	
    forAll(Sy_, i)
    {
        Sy_[i] = chemistry.RR(i);
    }

	NG_ = thermo_.p()/plasmaConstants::boltzC/thermo_.T();

	NG_.correctBoundaryConditions();

    return correct(chemistry, E, fields);
}

Foam::scalar 
Foam::multiSpeciesPlasmaModel::divFe()
{

	tmp<volScalarField> tdivFe
    (
        new volScalarField
        (
            IOobject
            (
                "tdivFe",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );


    volScalarField& divFe = tdivFe();

    //Info << "div step " << endl;

    divFe = mag(0.5*fvc::div(F_[eIndex_]));
    tdivFe().correctBoundaryConditions();

    //Info << "div step done " << endl;

    scalar maxdivFe = gMax(divFe);

    return maxdivFe;

}

inline Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::RR
(
    const word Yname,
    const psiChemistryModel& chemistry,
    const label i
) const
{
    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + Yname + ')',
                chemistry.time().timeName(),
                chemistry.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            chemistry.mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

	scalarField& Su = tRR();

    if (chemistry.chemistry())
    {
        forAll(Su, celli)
        {
            Su[celli] = chemistry.RR(i)()[celli];
        }

        tRR().correctBoundaryConditions();
    }
    return tRR;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::electronTempSource
(
	const psiChemistryModel& chemistry
)
{
    tmp<volScalarField> tEts
    (
        new volScalarField
        (
            IOobject
            (
                "electronTempSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    volScalarField& ets = tEts();

    volScalarField pets = 0.0*ets;

    pets = (chemistry.eChemSource()())*plasmaConstants::A/W(eIndex_); 

	if (collisionFrequency_[eIndex_] == "muBased")
	{
		pets += (3*plasmaConstants::boltzC*plasmaConstants::eChargeA*(thermo_.Te()-thermo_.T())*N(eIndex_)/mu_[eIndex_]/W(bIndex_));
	}

    forAll(ets, celli)
    {
        ets[celli] = pets[celli];
    }

    forAll(ets.boundaryField(), patchi)
    {
        const fvPatchScalarField& ppets =
            pets.boundaryField()[patchi];

        fvPatchScalarField& pppets = ets.boundaryField()[patchi];

        forAll(pppets, facei)
        {
            pppets[facei] = ppets[facei];
        }
    }
   
    return tEts;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::dElectronTempSourceDTe
(
	const psiChemistryModel& chemistry
)
{
    tmp<volScalarField> tEts
    (
        new volScalarField
        (
            IOobject
            (
                "dElectronTempSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    volScalarField& ets = tEts();

    volScalarField pets = 0.0*ets;

    pets = (chemistry.dEChemSourceDTe()())*plasmaConstants::A/W(eIndex_); 

	if (collisionFrequency_[eIndex_] == "muBased")
	{
		pets += (3*plasmaConstants::boltzC*plasmaConstants::eChargeA*N(eIndex_)/mu_[eIndex_]/W(bIndex_));
	}

    forAll(ets, celli)
    {
        ets[celli] = pets[celli];
    }

    forAll(ets.boundaryField(), patchi)
    {
        const fvPatchScalarField& ppets =
            pets.boundaryField()[patchi];

        fvPatchScalarField& pppets = ets.boundaryField()[patchi];

        forAll(pppets, facei)
        {
            pppets[facei] = ppets[facei];
        }
    }
   
    return tEts;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::ionTempSource
(
	const psiChemistryModel& chemistry,
	const volVectorField& E
)
{
    tmp<volScalarField> tIts
    (
        new volScalarField
        (
            IOobject
            (
                "ionTempSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    volScalarField& its = tIts();

    forAll(species(), i)
    {
        if (i != eIndex_ && i < activeSpecies_)
        {
            its += plasmaConstants::eCharge*(J(i) & E);

            if (collisionFrequency_[i] == "muBased")
            {
                its += (0.75*plasmaConstants::boltzC*plasmaConstants::eChargeA*(thermo_.Tion()-thermo_.T())*N(i)/mu_[i]/W(i));
            }
            else if (collisionFrequency_[i] == "rrBased")
            {
                its += (0.75*plasmaConstants::boltzC*(thermo_.Tion()-thermo_.T())*N(i)*chemistry.collFreq(i));
            }
		}        
    }

    tIts().correctBoundaryConditions();
   
    return tIts;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::netCharge()
{
    tmp<volScalarField> tPs
    (
        new volScalarField
        (
            IOobject
            (
                "netCharge",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    volScalarField& ps = tPs();

    volScalarField pps = 0.0*ps;
    
    forAll(species(), i)
    {
        if (i < activeSpecies_)
        {	
            pps += z_[i]*N(i);
		}        
    }    

    forAll(ps, celli)
    {
        ps[celli] = pps[celli];
    }

    forAll(ps.boundaryField(), patchi)
    {
        const fvPatchScalarField& ppps =
            pps.boundaryField()[patchi];

        fvPatchScalarField& pppps = ps.boundaryField()[patchi];

        forAll(pppps, facei)
        {
            pppps[facei] = ppps[facei];
        }
    }

    return tPs;
}


Foam::tmp<Foam::volVectorField>
Foam::multiSpeciesPlasmaModel::gradTe()
{
    tmp<volVectorField> tgradTe
    (
        new volVectorField
        (
            IOobject
            (
                "gradTe",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
        )
    );

    volVectorField& gradTe = tgradTe();

    gradTe = fvc::grad((thermo_.Te()));

    return tgradTe;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::potentialImpSource()
{
    tmp<volScalarField> tPis
    (
        new volScalarField
        (
            IOobject
            (
                "potentialImpSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    volScalarField& pis = tPis();
    
    forAll(species(), i)
    {
        if (i < activeSpecies_)
        {
		    volScalarField netImpSource = sign(z_[i])*z_[i]*mu_[i]*N(i);
		
		    forAll(pis, celli)
		    {
		        pis[celli] += netImpSource[celli];
		    }

		    forAll(pis.boundaryField(), patchi)
		    {
		        const fvPatchScalarField& pnetImpSource =
		            netImpSource.boundaryField()[patchi];

		        fvPatchScalarField& ppis = pis.boundaryField()[patchi];

		        forAll(ppis, facei)
		        {
		            ppis[facei] += pnetImpSource[facei];
		        }
		    }
		}        
    }    
    return tPis;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::potentialExpSource()
{
    tmp<volScalarField> tPes
    (
        new volScalarField
        (
            IOobject
            (
                "potentialExpSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    volScalarField& pes = tPes();
    
    forAll(species(), i)
    {
        if (i < activeSpecies_)
        {
        	//Info << "N(i) old = " << N(i).oldTime() << endl;
		    volScalarField netExpSource 
			= z_[i]*(N(i).oldTime()
			+ runTime_.deltaT().value()*((Sy_[i]*plasmaConstants::A/W(i))
			+ fvc::laplacian(D_[i], N(i), "laplacian(D,Ni)")
			+ fvc::laplacian((mu_[i]*plasmaConstants::KBE*N(i)), T_[i], "laplacian(D,T)")));
		
		    forAll(pes, celli)
		    {
		        pes[celli] += netExpSource[celli];
		    }
		}        
    }    
    return tPes;
}

Foam::tmp<Foam::volVectorField>
Foam::multiSpeciesPlasmaModel::netChargeFlux()
{
    tmp<volVectorField> tNcf
    (
        new volVectorField
        (
            IOobject
            (
                "netChargeFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(1, -1, -3, 0, 0), vector::zero)
        )
    );
    
    volVectorField& ncf = tNcf();

	volVectorField pNcf = 0.0*ncf;
    
    forAll(species(), i)
    {
        if (i < activeSpecies_)
        {
            pNcf += z_[i]*J(i);
		}        
    }

    pNcf*= plasmaConstants::eCharge;

    forAll(ncf, celli)
    {
        ncf[celli] = pNcf[celli];
    }

    forAll(ncf.boundaryField(), patchi)
    {
        const fvPatchVectorField& ppncf =
            pNcf.boundaryField()[patchi];

        fvPatchVectorField& pppncf = ncf.boundaryField()[patchi];

        forAll(pppncf, facei)
        {
            pppncf[facei] = ppncf[facei];
        }
    }

    return tNcf;
}



Foam::tmp<Foam::volVectorField>
Foam::multiSpeciesPlasmaModel::totalIonFlux()
{
    tmp<volVectorField> tIf
    (
        new volVectorField
        (
            IOobject
            (
                "totalIonFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(1, -1, -3, 0, 0), vector::zero)
        )
    );
    
    volVectorField& If = tIf();

	volVectorField pIf = 0.0*If;
    
    forAll(species(), i)
    {
        if (i != eIndex_ && i < activeSpecies_)
        {
            pIf += J(i);
		}        
    }

    forAll(If, celli)
    {
        If[celli] = pIf[celli];
    }

    forAll(If.boundaryField(), patchi)
    {
        const fvPatchVectorField& ppIf =
            pIf.boundaryField()[patchi];

        fvPatchVectorField& pppIf = If.boundaryField()[patchi];

        forAll(pppIf, facei)
        {
            pppIf[facei] = ppIf[facei];
        }
    }

    return tIf;
}

Foam::tmp<Foam::volVectorField>
Foam::multiSpeciesPlasmaModel::electronConvectiveFlux()
{
    tmp<volVectorField> tEcf
    (
        new volVectorField
        (
            IOobject
            (
                "electronConvectiveFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(1, -1, -3, 0, 0), vector::zero)
        )
    );
    
    volVectorField& ecf = tEcf();
    
    volVectorField pEcf = F_[eIndex_];

    forAll(ecf, celli)
    {
        ecf[celli] = pEcf[celli];
    }

    forAll(ecf.boundaryField(), patchi)
    {
        const fvPatchVectorField& ppEcf =
            pEcf.boundaryField()[patchi];

        fvPatchVectorField& pppEcf = ecf.boundaryField()[patchi];

        forAll(pppEcf, facei)
        {
            pppEcf[facei] = ppEcf[facei];
        }
    }

    return tEcf;
}

Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesPlasmaModel::electronConductivity
(
	const psiChemistryModel& chemistry
)
{
    tmp<volScalarField> tEc
    (
        new volScalarField
        (
            IOobject
            (
                "electronConductivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    volScalarField& ec = tEc();

	volScalarField collFreqT = 0.0*ec;

	if (collisionFrequency_[eIndex_] == "rrBased")
	{
		collFreqT = chemistry.collFreq(eIndex_)*W(eIndex_)*plasmaConstants::rA;
    }
	else if (collisionFrequency_[eIndex_] == "muBased")
	{
		collFreqT = plasmaConstants::eCharge/mu_[eIndex_];
	}

	volScalarField ptEc = 2.5*plasmaConstants::boltzCsqr*T_[eIndex_]*N(eIndex_)/collFreqT;

    forAll(ec, celli)
    {
        ec[celli] = ptEc[celli];
    }

    forAll(ec.boundaryField(), patchi)
    {
        const fvPatchScalarField& pptEc =
            ptEc.boundaryField()[patchi];

        fvPatchScalarField& ppEc = ec.boundaryField()[patchi];

        forAll(ppEc, facei)
        {
            ppEc[facei] = pptEc[facei];
        }
    }
    return tEc;
}

bool Foam::multiSpeciesPlasmaModel::read()
{
    return regIOobject::read();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
