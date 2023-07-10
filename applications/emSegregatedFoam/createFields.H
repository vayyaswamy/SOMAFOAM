// Block vector field for magnetic vector potential real part (first entry) and imaginary part (second entry)
    Info << "Creating field A\n" << endl;
    volVector6Field A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector6("zero", dimless, vector6::zero)
    );

    Info<< "Reading dictionary: physicalProperties\n" << endl;

// read .../constant/physicalerProperties
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar muMag
    (
        physicalProperties.lookup("muMag")
    );


    dimensionedScalar w
    (
        physicalProperties.lookup("w")
    );

    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    
    volScalarField DA
    (
        IOobject
        (
            "DA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("DA",dimensionSet(0,-2,0,0,0,0,0),0)
    );

    DA = muMag*sigma*w; // Auxiliary variable 

    Info<< "Reading intial conditions of field variables from .../0 directory\n" << endl;

    volVectorField Jcoil
    (
        IOobject
        (
            "Jcoil",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField AR
    (
        IOobject
        (
        "AR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField deltaAR
    (
        IOobject
        (
        "deltaAR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,1,-2,0,0,-1,0), vector::zero)
    );

    // CQ: AR_last to save AR from last iteration, to be used in relaxation
    volVectorField AR_last
    (
        IOobject
        (
        "AR_last",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        AR                 // CQ: initialization, not important
    );


    volVectorField AI
    (
        IOobject
        (
        "AI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField deltaAI
    (
        IOobject
        (
        "deltaAI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,1,-2,0,0,-1,0), vector::zero)
    );

    // CQ: AI_last to save AR from last iteration, to be used in relaxation
    volVectorField AI_last
    (
        IOobject
        (
        "AI_last",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        AI                 // CQ: initialization not important
    );

    volScalarField Amag
    (
        IOobject
        (
            "Amag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((AR & AR) + (AI & AI)) // A = sqrt(AR*AR + AI*AI) 
    );

    volVectorField BR
    (
        IOobject
        (
            "BR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(AR) 
    );

    volVectorField BI
    (
        IOobject
        (
            "BI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(AI)
    );

    volScalarField Bmag
    (
        IOobject
        (
            "Bmag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((BR & BR) + (BI & BI))
    );

    volScalarField Hmag
    (
        IOobject
        (
            "Hmag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Bmag/muMag     
    );

    volVectorField JR
        (
            IOobject
            (
            "JR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sigma*w*AI
    );

    volVectorField JI
        (
            IOobject
            (
            "JI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -sigma*w*AR
    );

    volScalarField Jind
    (
        IOobject
        (
			"Jind",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
        ),
        sqrt((JR & JR) + (JI & JI))   
    );

	volScalarField qJ
	(
		IOobject
		(
			"qJ",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		sigma/2.0 * sqr(w) * sqr(Amag)  
	);

    volVectorField fL
        (
            IOobject
            (
            "fL",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5 * ((JR ^ BR) - (JI ^ BI))
    );

