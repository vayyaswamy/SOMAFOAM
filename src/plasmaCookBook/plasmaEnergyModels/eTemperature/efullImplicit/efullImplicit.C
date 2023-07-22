/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "efullImplicit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(efullImplicit, 0);
    addToRunTimeSelectionTable(eTemp, efullImplicit, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efullImplicit::efullImplicit
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E,
    const dictionary& dict
)
:
    eTemp(thermo, mspm, E),
    eeFlux
    (
        IOobject
        (
            "eeFlux",
            thermo.T().mesh().time().timeName(),
			thermo.T().mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        thermo.T().mesh(),
        dimensionedVector("zero", dimensionSet(0, 0, 0, 1, 0), vector::zero)
    ),
    eSpecie("electron"),
    eIndex_(mspm.species()[eSpecie])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::efullImplicit::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E
)
{

    lduSolverPerformance solverPerf;
	volScalarField& TeC = thermo().Te();
    TeC.storePrevIter();
    scalar initialResidual = 1.0;

    int icorr = 0;

    while ((icorr++ <= 5) & (initialResidual >= 1e-5))

    {

	   eeFlux = 2.5*plasmaConstants::boltzC*mspm().J(eIndex_);

	   surfaceScalarField eeFluxF = fvc::interpolate(eeFlux) & mesh().Sf();

	   volScalarField eeSource = - plasmaConstants::eCharge*(mspm().J(eIndex_) & E) - mspm().electronTempSource(chemistry);

        volScalarField eeSource_Su = plasmaConstants::eCharge*(mspm().J(eIndex_) & E) 
                                    + mspm().electronTempSource(chemistry) 
                                    - mspm().dElectronTempSourceDTe(chemistry)*TeC;

        volScalarField eeSource_SuSp = mspm().dElectronTempSourceDTe(chemistry);

        //Info << eeSource_SuSp << endl;

	   const volScalarField& Ne = mspm().N(eIndex_);

        fvScalarMatrix TeEqn
        (
            fvm::ddt((1.5*plasmaConstants::boltzC*Ne), TeC)
            + fvm::div(eeFluxF, TeC)
            - fvm::laplacian(mspm().electronConductivity(chemistry), TeC, "laplacian(eC,Te)")
            + fvm::SuSp(eeSource_SuSp, TeC)
            + eeSource_Su 
        );

	   solverPerf = TeEqn.solve();

        initialResidual = solverPerf.initialResidual();

        TeC.max(300.0);

        thermo().correct();

        mspm().updateTemperature();

        mspm().updateChemistryCollFreq(chemistry);

    }

    TeC.relax();

    TeC.max(300.0);

    thermo().correct();

    mspm().updateTemperature();

    mspm().updateChemistryCollFreq(chemistry);

    volScalarField divEeFlux = mag(0.5*fvc::div(eeFlux));

    scalar maxDivEeFlux = gMax(divEeFlux);

    return maxDivEeFlux;
}


// ************************************************************************* //
