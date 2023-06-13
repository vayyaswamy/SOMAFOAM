/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "fullImplicitTE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(efullImplicitTE, 0);
    addToRunTimeSelectionTable(eTemp, efullImplicitTE, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efullImplicitTE::efullImplicitTE
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
			IOobject::MUST_READ,
			IOobject::NO_WRITE
        ),
        thermo.T().mesh()
    ),
    eSpecie("electron"),
    eIndex_(mspm.species()[eSpecie]),
    fieldBounds(thermo.T().mesh().solutionDict().subDict("fieldBounds")),
    TeMin(0.0),
    TeMax(0.0)
{
    fieldBounds.lookup(thermo.Te().name()) >> TeMin >> TeMax;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::efullImplicitTE::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E
)
{
	volScalarField& TeC = thermo().Te();

	eeFlux = 2.5*plasmaConstants::boltzC*mspm().J(eIndex_)*TeC;

	eeFlux.correctBoundaryConditions();

	if (restartCapable() && runTime().write())
	{
		eeFlux.write();
	}

	surfaceScalarField eeFluxF = fvc::interpolate(eeFlux/TeC) & mesh().Sf();

	surfaceScalarField eFluxF = fvc::interpolate(mspm().J(eIndex_)) & mesh().Sf();

	volScalarField eeSource = -plasmaConstants::eCharge*(mspm().J(eIndex_) & E) - mspm().electronTempSource(chemistry);

	const volScalarField& Ne = mspm().N(eIndex_);

    volScalarField K = 0.5*plasmaConstants::rA*mspm().W(eIndex_)*magSqr(mspm().U(eIndex_));

    fvScalarMatrix TeEqn
    (
        fvm::ddt((1.5*plasmaConstants::boltzC*Ne), TeC)
      + fvm::div(eeFluxF, TeC, "div(eeFlux,Te)")
      + fvc::ddt(Ne, K)
      + fvc::div(eFluxF, K, "div(eFlux,K)") 
      - fvm::laplacian(mspm().electronConductivity(chemistry), TeC, "laplacian(eC,Te)")
	  + fvm::SuSp((-eeSource/TeC), TeC)
    );

    TeEqn.relax();

	TeEqn.solve();

    TeC.max(TeMin);

	TeC.min(TeMax);

    return 0;
}


// ************************************************************************* //
