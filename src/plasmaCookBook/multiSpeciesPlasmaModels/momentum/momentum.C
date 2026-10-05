/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "momentum.H"
#include "zeroGradientFvPatchFields.H"
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
inline void Foam::momentum<ThermoType>::updateVelocity
(
    const psiChemistryModel& chemistry,
    const volVectorField& E,
    const label i
)
{
    const volScalarField& Ti = T_[i];

    const volScalarField& Ni = N_[i];

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

    // Momentum equation divided by the particle mass m = W/A
    tmp<fvVectorMatrix> uEqn
    (
        fvm::ddt(Ni, Ui)
      + pConvection->fvmDiv(Jif, Ui)
      + fvm::Sp((Ni*collFrequency_[i]), Ui)
     ==
        (
          - fvc::grad(plasmaConstants::boltzC*Ti*Ni)
          + z_[i]*plasmaConstants::eCharge*Ni*E
        )*plasmaConstants::A/W(i)
    );

    uEqn->relax();

    uEqn->solve(mesh_.solutionDict().solver("Ui"));

    if (restartcapabale && runTime_.outputTime())
    {
        Ui.write();
    }
}

template<class ThermoType>
inline void Foam::momentum<ThermoType>::transportCoeffInterpolate
(
    const label i
)
{
    if (i < activeSpecies_ && collisionFrequency_[i] == "muBased")
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
    else
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
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::momentum<ThermoType>::momentum
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
    U_.setSize(activeSpecies_);

    // F_ mirrors U_ so that divFe() and electronConvectiveFlux() work
    F_.setSize(activeSpecies_);

    updateTemperature();

    forAll(U_, i)
    {
        const word Uname("U_" + species()[i]);

        IOobject header
        (
            Uname,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        );

        if (header.headerOk())
        {
            U_.set
            (
                i, new volVectorField
                (
                    IOobject
                    (
                        Uname,
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
            U_.set
            (
                i, new volVectorField
                (
                    IOobject
                    (
                        Uname,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("zero", dimVelocity, vector::zero),
                    zeroGradientFvPatchVectorField::typeName
                )
            );
        }

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
                U_[i]
            )
        );
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::momentum<ThermoType>::solveSpecie
(
    const label i,
    psiChemistryModel& chemistry,
    const volVectorField& E
)
{
    volScalarField& yi = thermo_.composition().Y(i);
    volScalarField& Ni = N_[i];

    transportCoeffInterpolate(i);

    if (i < activeSpecies_)
    {
        updateVelocity(chemistry, E, i);

        fvScalarMatrix NEqn
        (
            fvm::ddt(Ni)
          + fvm::div((fvc::interpolate(U_[i]) & mesh_.Sf()), Ni, "div(U,Ni)")
          + fvm::SuSp((-Sy_[i]*plasmaConstants::A/W(i)/Ni), Ni)
        );

        NEqn.relax();

        NEqn.solve(mesh_.solutionDict().solver("Ni"));

        this->limitField(Ni, 1e6);

        J_[i] == Ni*U_[i];

        F_[i] == U_[i];
    }
    else
    {
        fvScalarMatrix NnEqn
        (
            fvm::ddt(Ni)
          - fvm::laplacian(D_[i], Ni, "laplacian(D,Nin)")
          + fvm::SuSp((-Sy_[i]*plasmaConstants::A/W(i)/Ni), Ni)
        );

        NnEqn.relax();

        NnEqn.solve(mesh_.solutionDict().solver("Nin"));

        this->limitField(Ni, 1e4);
    }

    yi = Ni*W(i)/thermo_.rho()/plasmaConstants::A;
}


template<class ThermoType>
inline Foam::scalar Foam::momentum<ThermoType>::correct
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
bool Foam::momentum<ThermoType>::read()
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
