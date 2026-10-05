/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "driftDiffusion.H"
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
inline void Foam::driftDiffusion<ThermoType>::updateFlux
(
    const label i,
    const volVectorField& E
)
{
    volVectorField& Fi = F_[i];

    // Drift velocity (electric field + thermal drift) per unit density
    Fi = mu_[i]*(z_[i]*E - plasmaConstants::KBE*fvc::grad(T_[i]));

    if (restartcapabale && runTime_.outputTime())
    {
        Fi.write();
    }
}

template<class ThermoType>
inline void Foam::driftDiffusion<ThermoType>::transportCoeffInterpolate
(
    const label i
)
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

    if ( i < activeSpecies_)
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
        if(diffusionModel_[i] == "einsteinRelation")
        {
            D_[i] = mu_[i]*plasmaConstants::KBE*T_[i];

            D_[i].correctBoundaryConditions();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::driftDiffusion<ThermoType>::driftDiffusion
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
    F_.setSize(activeSpecies_);

    updateTemperature();

    forAll(F_, i)
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
                dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector::zero)
            )
        );
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::driftDiffusion<ThermoType>::solveSpecie
(
    const label i,
    psiChemistryModel& chemistry,
    const volVectorField& E
)
{
    volScalarField& yi = thermo_.composition().Y(i);
    volScalarField& Ni = N_[i];

    transportCoeffInterpolate(i);

    const bool charged = (i < activeSpecies_);

    if (charged)
    {
        updateFlux(i, E);
    }

    const word solverName(charged ? "Ni" : "Nin");
    const label nCorr = charged ? nCorrCharged_ : nCorrNeutral_;

    Ni.storePrevIter();

    // Chemistry source linearised in Ni and updated each corrector
    scalar initialResidual = 1.0;
    label iCorr = 0;

    while ((initialResidual >= innerTolerance_) && (iCorr++ < nCorr))
    {
        fvScalarMatrix NEqn
        (
            fvm::ddt(Ni)
          - chemistry.RR(i)*plasmaConstants::A/W(i)
          + chemistry.dRRDi(i)*Ni
          - fvm::Sp(chemistry.dRRDi(i), Ni)
        );

        if (charged)
        {
            NEqn +=
                fvm::div((fvc::interpolate(F_[i]) & mesh_.Sf()), Ni, "div(F,Ni)")
              - fvm::laplacian(D_[i], Ni, "laplacian(D,Ni)");
        }
        else
        {
            NEqn -= fvm::laplacian(D_[i], Ni, "laplacian(D,Nin)");
        }

        initialResidual =
            NEqn.solve(mesh_.solutionDict().solver(solverName)).initialResidual();

        if (!charged)
        {
            this->limitField(Ni, densityFloor_);
        }

        yi = Ni*W(i)/thermo_.rho()/plasmaConstants::A;

        updateChemistryCollFreq(chemistry);
    }

    reportInnerIterations(species()[i], nCorr, initialResidual);

    relaxFinal(Ni);

    this->limitField(Ni, densityFloor_);

    yi = Ni*W(i)/thermo_.rho()/plasmaConstants::A;

    updateChemistryCollFreq(chemistry);

    if (charged)
    {
        const volScalarField& Ti = T_[i];

        if (diffusionModel_[i] == "einsteinRelation")
        {
            J_[i] == -mu_[i]*
            (
                Ti*fvc::grad(plasmaConstants::KBE*Ni)
              + plasmaConstants::KBE*Ni*fvc::grad(Ti)
              - Ni*z_[i]*E
            );
        }
        else
        {
            J_[i] ==
              - D_[i]*fvc::grad(Ni)
              - plasmaConstants::KBE*mu_[i]*Ni*fvc::grad(Ti)
              + Ni*mu_[i]*z_[i]*E;
        }
    }
}


template<class ThermoType>
inline Foam::scalar Foam::driftDiffusion<ThermoType>::correct
(
    psiChemistryModel& chemistry,
    const volVectorField& E,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{
    updateTemperature();

    if (eonCalculation)
    {
        EON = mag(E)*1E21/NG_;

        EON.correctBoundaryConditions();
    }

    volScalarField yt = 0.0*thermo_.composition().Y(0);

    forAll(species(), i)
    {
        if (i != bIndex_ && speciesSolution_[i])
        {
            if (multiTimeStep)
            {
                if ((timeS_[i] == "LTS") && (fmod(runTime_.deltaT().value(),LTScounter) < SMALL))
                {
                    LTSset();
                }
                else if ((timeS_[i] == "MTS") && (fmod(runTime_.deltaT().value(),MTScounter) < SMALL))
                {
                    MTSset();
                }
            }

            solveSpecie(i, chemistry, E);

            yt += thermo_.composition().Y(i);

            if (multiTimeStep)
            {
                STSset();
            }
        }
    }

    volScalarField& yBgas = thermo_.composition().Y()[bIndex_];

    yBgas = scalar(1.0) - yt;

    N_[bIndex_] == thermo_.rho()*yBgas*plasmaConstants::A/W(bIndex_);

    updateChemistryCollFreq(chemistry);

    return 0;
}

template<class ThermoType>
bool Foam::driftDiffusion<ThermoType>::read()
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
