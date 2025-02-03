/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "gfullImplicit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gfullImplicit, 0);
    addToRunTimeSelectionTable(gTemp, gfullImplicit, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gfullImplicit::gfullImplicit
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E,
    const dictionary& dict
)
:
    gTemp(thermo, mspm, E),
    tavesource(dict.lookup("mode")),
    cycleavevalue(dict.lookupOrDefault("cycleavevalue", 0.0)),
    gasTempSource
    (
        IOobject
        (
            "gasTempSource",
            thermo.T().mesh().time().timeName(),
			thermo.T().mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        thermo.T().mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
    ),
    gasTempSourcetemp
    (
        IOobject
        (
            "gasTempSourcetemp",
            thermo.T().mesh().time().timeName(),
			thermo.T().mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        thermo.T().mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
    ),
    timeCount(0.0),
    timeCountn(0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gfullImplicit::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E
)
{
	volScalarField& Tc = thermo().T();

    //if(tavesource == "frequency")
    //{
    //    gasTempSourcetemp += mspm().ionTempSource(chemistry, E) + chemistry.Sh()();
    //    timeCountn++;
    //    if(runTime().value()-timeCount > cycleavevalue)
		//{
    //        timeCount = runTime().value();
    //        gasTempSource = gasTempSourcetemp/timeCountn;
    //        gasTempSourcetemp = 0.0*gasTempSource;
    //        timeCountn = 0.0;
    //    }
    //}

    gasTempSource = mspm().ionTempSource(chemistry, E);
    volScalarField kappa = thermo().Cp()*thermo().rho()*thermo().alpha();
    volScalarField rhoCp = thermo().Cp()*thermo().rho();
    
    //Info << "Cp = " << thermo().Cp() << endl;
    //Info << "rho = " << thermo().rho() << endl;
    //Info << "alpha = " << thermo().alpha() << endl;
    //Info << "chemistry = " << chemistry.Sh() << endl;


    fvScalarMatrix TEqn
    (
      - fvm::laplacian(kappa, Tc, "laplacian(alpha,T)")
	  ==
      gasTempSource
    );

    TEqn.relax();

	TEqn.solve();
 
    Tc.relax();
}


// ************************************************************************* //
