    Info<< "Reading physicalProperties\n" << endl;

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

    dimensionedScalar me
    (
        physicalProperties.lookup("me")
    );
    dimensionedScalar m1
    (
        physicalProperties.lookup("m1")
    );
    dimensionedScalar m2
    (
        physicalProperties.lookup("m2")
    );

    dimensionedScalar Te
    (
        physicalProperties.lookup("Te")
    );

    dimensionedScalar T1
    (
        physicalProperties.lookup("T1")
    );

    dimensionedScalar T2
    (
        physicalProperties.lookup("T2")
    );
    
    dimensionedScalar Se
    (
        physicalProperties.lookup("Se")
    );

    dimensionedScalar S1
    (
        physicalProperties.lookup("S1")
    );

    dimensionedScalar S2
    (
        physicalProperties.lookup("S2")
    );

    Info<< "Reading field phi\n" << endl;
    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field ne\n" << endl;
    volScalarField ne
    (
        IOobject
        (
            "ne",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field n1\n" << endl;
    volScalarField n1
    (
        IOobject
        (
            "n1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    
    Info<< "Reading field n2\n" << endl;
    volScalarField n2
    (
        IOobject
        (
            "n2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field ni\n" << endl;
    volScalarField ni
    (
        IOobject
        (
            "ni",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        n1+n2
    );


    Info<< "Reading field mome\n" << endl;
    volVectorField Ue
    (
        IOobject
        (
            "Ue",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field mom1\n" << endl;
    volVectorField U1
    (
        IOobject
        (
    	    "U1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );	

    Info<< "Reading field mom2\n" << endl;
    volVectorField U2
    (
        IOobject
        (
    	    "U2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );	

    Info<< "Calculating field E\n" << endl;
    volVectorField E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -fvc::grad(phi)
    );


    Info<< "Calculating field Gammae\n" << endl;
    volVectorField Gammae
    (
        IOobject
	(
	    "Gammae",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE 
	),
	ne*Ue

    );

    Info<<"Calculating field Gamma1\n" << endl;
    volVectorField Gamma1
    (
        IOobject
	(
	     "Gamma1",
	     runTime.timeName(),
	     mesh,
	     IOobject::NO_READ,
	     IOobject::AUTO_WRITE
	),
	n1*U1
    );

    Info<<"Calculating field Gamma2\n" << endl;
    volVectorField Gamma2
    (
        IOobject
	(
	     "Gamma2",
	     runTime.timeName(),
	     mesh,
	     IOobject::NO_READ,
	     IOobject::AUTO_WRITE
	),
	n2*U2
    );
    Info<< "Calculating field neFlux\n" << endl;
    surfaceScalarField neFlux
    (
        IOobject
	(
	    "neFlux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	linearInterpolate(Ue) & mesh.Sf()
    );

    Info<< "Calculating field n1Flux\n" << endl;
    surfaceScalarField n1Flux
    (
        IOobject
	(
	    "n1Flux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	linearInterpolate(U1) & mesh.Sf()
    );

    Info<< "Calculating field n2Flux\n" << endl;
    surfaceScalarField n2Flux
    (
        IOobject
	(
	    "n2Flux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	linearInterpolate(U2) & mesh.Sf()
    );


    surfaceScalarField UeFlux
    (
        IOobject
	(
	    "UeFlux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	linearInterpolate(ne)*neFlux
    );

    Info<< "Calculating field U1Flux\n" << endl;
    surfaceScalarField U1Flux
    (
        IOobject
	(
	    "U1Flux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	linearInterpolate(n1)*n1Flux
    );

    Info<< "Calculating field U2Flux\n" << endl;
    surfaceScalarField U2Flux
    (
        IOobject
	(
	    "U2Flux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	linearInterpolate(n2)*n2Flux
    );
