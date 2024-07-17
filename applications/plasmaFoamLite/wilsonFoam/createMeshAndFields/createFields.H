//==============================================================================================
//----------------- READING FROM CONSTANT/REGIONPROPERTIES FILE  -------------------------------
//----------------------------------------------------------------------------------------------
IOdictionary dielectricRegions
    (
        IOobject
        (
            "regionProperties",
            runTime.time().constant(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
	//scalar convergenceDelta(readScalar(dielectricRegions.lookup("poissonConvergenceDeltaValue")));
	//scalar maxIter(readScalar(dielectricRegions.lookup("maxIterations")));
	//scalar secondaryElectronEmissionCoefficient(readScalar(dielectricRegions.lookup("secondaryElectronEmissionCoefficient")));
	Foam::List<Foam::word> dielectricList(dielectricRegions.lookup("dielectricRegionNames"));
	//Foam::List<Foam::word> inletAndOutletBoundary(dielectricRegions.lookup("inletAndOutletBoundaryNames"));
	//Foam::List<Foam::word> metallicBoundary(dielectricRegions.lookup("metallicBoundaryNames"));
	////Info << "Number of dielectric regions mentioned ::::  " << dielectricList.size()<< endl << endl;
	//Info << "Number of inlet and outlet Boundaries ::::  " << inletAndOutletBoundary.size()<< endl << endl;
	//Info << "Number of Metallic Boundaries ::::  " << metallicBoundary.size()<< endl << endl;
	//Info<<  "poissonConvergenceDeltaValue :::: " << convergenceDelta << endl;
	//Info<<  "secondaryElectronEmissionCoefficient :::: " << secondaryElectronEmissionCoefficient << endl;


scalar maxCo(readScalar(runTime.controlDict().lookup("maxCo")));
scalar deltaTMax(readScalar(runTime.controlDict().lookup("deltaTMax")));
scalar deltaTMin(readScalar(runTime.controlDict().lookup("deltaTMin")));
//scalar waitUntil(readScalar(runTime.controlDict().lookup("waitUntil")));

// =========================================================================
// ============= CREATE MESH ===============================================
Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;
// REf: Cht MultiRegionFoam tutorial. 
fvMesh mesh
(
    IOobject
    (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
//===========================================================================
PtrList<fvMesh> dielectricMesh(dielectricList.size());
forAll(dielectricList, i)
{
    dielectricMesh.set
    (
        i,
        new fvMesh
        (
            IOobject
            (
                dielectricList[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
//============================================================================
//============================== FIELDS ======================================
//============================================================================
const dimensionedScalar eCharge  //--------- ELIMENTARY CHARGE ---------------
(
    "eCharge",
    dimensionSet(0, 0, 1, 0, 0, 1, 0),
    1.6E-19
);
//========================================================================
const dimensionedScalar kBoltz  //--------- BOLTZMANN CONSTANT -----------
(
    "kBoltz",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    1.38065E-23
);
//========================================================================
 //-- THIS ARE USED TO CORRECT THE DIMENSION WHILE CALCULATING SOURCE TERM-
const dimensionedScalar dimensionCorrection
(
    "dimensionCorrection",
    dimensionSet(0, -3, 0, 0, 0, 0, 0),
    1
);

const dimensionedScalar dimensionCorrectionPower
(
    "dimensionCorrectionPower",
    dimensionSet(1, 2, -3, 0, 0, 0, 0),
    1
);
const dimensionedScalar dimensionCorrectionEPotential
(
    "dimensionCorrectionEPotential",
    dimensionSet(1, 2, -3, 0, 0, -1, 0),
    1
);

const dimensionedScalar dimensionCorrectionMobility
(
    "dimensionCorrectionMobility",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    1
);

const dimensionedScalar dimensionCorrectionEnergy
(
    "dimensionCorrectionEnergy",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    1
);
const dimensionedScalar dimensionCorrectionTemp
(
    "dimensionCorrectionTemp",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    1
);

const dimensionedScalar dimensionCorrectionTime
(
    "dimensionCorrectionTime",
    dimensionSet(0, 0, 1, 0, 0, 0, 0),
    1
);
const dimensionedScalar dimensionCorrectionEfield
(
    "dimensionCorrectionEfield",
    dimensionSet(1, 1, -3, 0, 0, -1, 0),
    1
);

const dimensionedScalar dimensionCorrectionMass
(
    "dimensionCorrectionMass",
    dimensionSet(1, 0, 0, 0, 0, 0, 0),
    1
);


//===============================================================================================
// ------------------------ NUMBER DENSITY OF ALL SPECIES ---------------------------------------
// ----------------------------------------------------------------------------------------------
PtrList<volScalarField> nDensity;
nDensity.setSize(speciesList.size());
forAll(speciesList, ii)
    {
		nDensity.set
        (
            ii, new volScalarField
            (
                IOobject
                (
                    speciesList[ii],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
	}
volScalarField unitMatrixx
(
	IOobject
	(
		"unitMatrixx",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1),
	zeroGradientFvPatchScalarField::typeName
);



//===============================================================================================
// ------------------------------- ELECTRON TEMPERATURE -----------------------------------------
// ----------------------------------------------------------------------------------------------
volScalarField eTemp
(
	IOobject
	(
		"eTemp",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
volScalarField eDensity
(
	IOobject
	(
		"eDensity",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh//,
	//dimensionedScalar("one", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0),
	//zeroGradientFvPatchScalarField::typeName
);
eDensity= (3/2)*kBoltz*eTemp*nDensity[0];
//eDensity.correctBoundaryConditions();
//-------------------------------------------------------------------------------------------------
volVectorField gradTe
(
		IOobject
		(
			"gradTe",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, -1, 0, 1, 0, 0, 0), vector::zero)
);
volVectorField gradne
(
		IOobject
		(
			"gradne",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, -4, 0, 0, 0, 0, 0), vector::zero)
);
volVectorField gradeD //gradient of energy density
(
		IOobject
		(
			"gradeD",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
);
//electron conductivity ------------------------------------------------------------
       volScalarField kappa_e
        (
                IOobject
                (
                    "kappa_electron",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimensionSet(1, 1,-3, -1, 0, 0, 0), 0.0),
				zeroGradientFvPatchScalarField::typeName
        );
// -----------------------------------------------------------------------------

// ============================ PRODUCTION AND SOURCE ==========================
volScalarField eSource
(
	IOobject
	(
		"eSource",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)//,
	//zeroGradientFvPatchScalarField::typeName
);
volScalarField totEnergyGen= 0*eSource;

//================================ JOULE HEATING ================================
volScalarField jouleHeating
(
	IOobject
	(
		"jouleHeating",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)//,
	//zeroGradientFvPatchScalarField::typeName
);


//==============================================================================================
// -------------------------------  TRANSPORT COEFFICIENTS  ------------------------------------
// ---------------------------------------------------------------------------------------------
PtrList<volScalarField> mobility;
PtrList<volScalarField> diffusion;
PtrList<volVectorField> velocity;
PtrList<volVectorField> Flux_;
PtrList<volScalarField> pSource;
PtrList<volScalarField> nMu;
PtrList<volScalarField> collFreqxN;

mobility.setSize(speciesList.size());
diffusion.setSize(speciesList.size());
velocity.setSize(speciesList.size());
Flux_.setSize(speciesList.size());
pSource.setSize(speciesList.size());
nMu.setSize(speciesList.size());
collFreqxN.setSize(speciesList.size());

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  MOBILITY  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                         
forAll (speciesList, ii){
	//collision frequency times ndensity
	collFreqxN.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"collFreqxN" + speciesList[ii],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)//,
			//zeroGradientFvPatchScalarField::typeName
		)
	);
	//----------------------------------------------------------
	mobility.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"mobility" + speciesList[ii],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), convertStringToFloat(speciesProperties[ii][4])/backgroundGasDensity),
			zeroGradientFvPatchScalarField::typeName
		)
	);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  DIFFUSION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	diffusion.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"diffusion" + speciesList[ii],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), convertStringToFloat(speciesProperties[ii][5])/backgroundGasDensity),
			zeroGradientFvPatchScalarField::typeName
		)
	);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  VELOCITY  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (transportModel[ii] =="momentum"){ 
		velocity.set
		(
			ii, new volVectorField
			(
				IOobject
				(
					"velocity" + speciesList[ii],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh//,
				//dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)//,
				//zeroGradientFvPatchVectorField::typeName
			)
		);
	}
	else{
		velocity.set
		(
			ii, new volVectorField
			(
				IOobject
				(
					"velocity" + speciesList[ii],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)//,
				//zeroGradientFvPatchVectorField::typeName
			)
		);
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FLUX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Flux_.set
		(
			ii, new volVectorField
			(
				IOobject
				(
					"Flux_" + speciesList[ii],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedVector("zero", dimensionSet(0, -2, -1, 0, 0, 0, 0), vector::zero)//,
				//zeroGradientFvPatchVectorField::typeName
			)
		);
	pSource.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"pSource" + speciesList[ii],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
	nMu.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"nMu" + speciesList[ii],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(-1, -3, 2, 0, 0, 1, 0), 0)
			//zeroGradientFvPatchScalarField::typeName
		)
	);
}


surfaceVectorField eFlux	
(
		IOobject
		(
			"eFlux",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0, 0, 0),  vector::zero)//,
		//zeroGradientFvPatchVectorField::typeName	
);



//==============================================================================================
// -------------------------------  REACTION COEFFICIENTS  -------------------------------------
// ---------------------------------------------------------------------------------------------
PtrList<volScalarField> reaction;
PtrList<volScalarField> reactionSource;
reaction.setSize(reactionCount);
reactionSource.setSize(reactionCount);

for (int ii=0; ii<reactionCount; ii++){
	reaction.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"reaction" + convertIntToString(ii),
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), reactionRate[ii])//,
			//zeroGradientFvPatchScalarField::typeName
		)
	);

	reactionSource.set
	(
		ii, new volScalarField
		(
			IOobject
			(
				"reactionSource" + convertIntToString(ii),
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)//,
			//zeroGradientFvPatchScalarField::typeName
		)
	);


}

//================================================================================================
//=================================== Energy Diffusion CHARGE ====================================
//================================================================================================

volScalarField eDiffusion
(
	IOobject
	(
		"eDiffusion",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);
//================================================================================================
//=================================== ELECTRON CONTUDCTIVITY =====================================
//================================================================================================

/*volScalarField kappa_e
(
	IOobject
	(
		"kappa_e",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);
*/
//================================================================================================
//======================================== SURFACE CHARGE ========================================
//================================================================================================

volScalarField surfC
(
	IOobject
	(
		"surfC",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh//,
	//dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), 0.0)//,
	//zeroGradientFvPatchScalarField::typeName
);

//================================================================================================
//========================================== NET CHARGE ==========================================
//================================================================================================

volScalarField netCharge
(
	IOobject
	(
		"netCharge",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, -3, 1, 0, 0, 1, 0), 0.0)//,
	//zeroGradientFvPatchScalarField::typeName
);
//================================================================================================
//======================= COLLISION FREQUENCY TIMES NUMBER DENSITY ===============================
//================================================================================================

volScalarField collisionFreqxN
(
	IOobject
	(
		"collisionFreqxN",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

//================================================================================================
//===================================== ELECTRIC POTENTIAL =======================================
//================================================================================================
volScalarField EPotential
(
	IOobject
	(
		"EPotential",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh//,
	//dimensionedScalar("zero", dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.0),
	//zeroGradientFvPatchScalarField::typeName

);
//================================================================================================
//========================================= POWER ================================================
//================================================================================================
volScalarField powerDissipated
(
	IOobject
	(
		"powerDissipated",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, 2, -3, 0, 0, 0, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);
volScalarField powerDissipatedElectron
(
	IOobject
	(
		"powerDissipatedElectron",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)//,
	//zeroGradientFvPatchScalarField::typeName
);
volScalarField powerDissipatedIon
(
	IOobject
	(
		"powerDissipatedIon",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)//,
	//zeroGradientFvPatchScalarField::typeName
);
volScalarField power //extra variable created for workaround for power BC
(
	IOobject
	(
		"power",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(1, 2, -3, 0, 0, 0, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

//================================================================================================
//===================================== DIELECTRIC CONSTANT=======================================
//================================================================================================
volScalarField epsilon
(
	IOobject
	(
		"epsilon",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh//,
	//dimensionedScalar("zero", dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.0),
	//zeroGradientFvPatchScalarField::typeName

);
//----- for calculation of impedance ---------------------------------------
/*volScalarField epsilon_impedance
(
	IOobject
	(
		"epsilon_impedance",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(-1, -3, 4, 0, 0, 2, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);
//--------------------------------------------------------------------------
volScalarField sigma_impedance
(
	IOobject
	(
		"epsilon_impedance",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(-1, -3, 5, 0, 0, 2, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);*/
volScalarField eCollisionFrequency
(
	IOobject
	(
		"eCollisionFrequency",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);
//================================================================================================
//================================ VELOCITY OF BACKGROUND GAS ====================================
//================================================================================================
volVectorField gasVelocity	
(
		IOobject
		(
			"gasVelocity",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh//,
		//dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(1, 0, 0)),
		//zeroGradientFvPatchVectorField::typeName
	
);
// ================== total ion flux ==============================================================
volVectorField ionFlux	
(
		IOobject
		(
			"ionFlux",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, -2, -1, 0, 0, 0, 0), vector(1, 0, 0))//,
		//zeroGradientFvPatchVectorField::typeName
	
);
// ================== total ion Density ==============================================================
volScalarField ionDensity	
(
		IOobject
		(
			"ionDensity",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0)//,
		//zeroGradientFvPatchVectorField::typeName
	
);

//================================================================================================
//================================ VELOCITY OF ELECTRON ENERGY ====================================
//================================================================================================
volVectorField eVelocity	
(
		IOobject
		(
			"eVelocity",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0, 0, 0))
		//zeroGradientFvPatchVectorField::typeName
	
);
surfaceVectorField tFlux= fvc::interpolate(2.5*kBoltz*eVelocity*nDensity[0]); //isnt exactly tFlux
//================================================================================================
//======================================== ELECTRIC FIELD ========================================
//================================================================================================
volVectorField EField
(
	IOobject
	(
		"EField",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)//,
	//zeroGradientFvPatchVectorField::typeName
);

//================================================================================================
//========================================== CELL VOLUME =========================================
//================================================================================================
volScalarField meshV
    (
        IOobject
        (
            "meshV",
            EField.mesh().time().timeName(),
            EField.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        EField.mesh(),
		dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0,0,0), 0.0)
    );
	meshV.internalField() = EField.mesh().V();

volScalarField EFieldTd = mag(EField)/backgroundGasDensity/dimensionCorrection/1E-21; 
volScalarField eFieldMag = mag(EField); 
surfaceScalarField avgNdensity=fvc::interpolate(nDensity[0]);
surfaceScalarField avgeTemp=fvc::interpolate(eTemp);
//=================================================================================================
volScalarField divFe=fvc::div(eFlux & mesh.Sf());
volScalarField potentialSource = chargeNumber[0]*eCharge*(fvc::div((fvc::interpolate(velocity[0]) & mesh.Sf()),nDensity[0],"div(vel,nDensity)") -
					fvc::laplacian(diffusion[0], nDensity[0],"laplacian(difC,nDensity)"));
volScalarField pTemperature= 0*eTemp;
//volVectorField pVel = 0*velocity[0];
volScalarField impSource = eCharge*chargeNumber[0]*nDensity[0]*mobility[0]*runTime.deltaT()/epsilon;
volVectorField totalElectronFlux= 0*eCharge*nDensity[0]*velocity[0]; //used for power calc
volVectorField totalIonFlux= 0*eCharge*nDensity[0]*velocity[0]; //used for power calc
//volVectorField totalIonFlux= nDensity[0]*velocity[0]; //used for secondary emission calc BC
volScalarField keterm1 = 0.5*nDensity[0]*(velocity[0] & velocity[0])*speciesMass[0]*dimensionCorrectionMass;
volVectorField keterm2 = 0.5*nDensity[0]*velocity[0]*(velocity[0]&velocity[0])*speciesMass[0]*dimensionCorrectionMass;

volVectorField currentDensity = Flux_[0]*eCharge;

/*surfaceScalarField phie
(
    IOobject
    (
        "phie",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (fvc::interpolate(velocity[0]) & mesh.Sf())
);*/
