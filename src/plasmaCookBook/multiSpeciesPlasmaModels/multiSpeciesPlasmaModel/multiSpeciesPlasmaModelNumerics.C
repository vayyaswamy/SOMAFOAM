/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Description
    Numerical controls of multiSpeciesPlasmaModel read from plasmaProperties:

    Inner iterations (defaults shown):
    \verbatim
    innerIterations
    {
        chargedSpecies      4;      // maximum passes per time step
        neutralSpecies      2;
        electronTemperature 6;
        tolerance           1e-5;   // initial residual that ends the passes
        reportUnconverged   yes;
    }
    \endverbatim

    Acceleration of slow neutral species such as metastables, whose
    diffusion and loss times are many periods of the applied voltage:
    \verbatim
    slowSpeciesAcceleration
    {
        species         (Arm);
        period          7.374631e-8;  // period of the applied voltage [s]
        fullCycles      5;      // periods of full simulation per block
        averageCycles   1;      // periods at the end of a block over which
                                // the sources are averaged
        deltaT          1e-6;   // time step of the advance [s]
        nSteps          100;    // steps per advance
        maxChangeFactor 2;      // limit on the density change per advance
        tolerance       1e-3;   // advances stop below this relative change
    }
    \endverbatim

    After each block the listed species alone are advanced by
    nSteps*deltaT with the period-averaged chemistry source, linearised in
    the species' own density, and diffusion; the charged species, the
    electron temperature and the potential are left as they are and adjust
    during the next block.

\*---------------------------------------------------------------------------*/

#include "multiSpeciesPlasmaModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fvm.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multiSpeciesPlasmaModel::readNumericalControls()
{
    if (found("innerIterations"))
    {
        const dictionary& dict = subDict("innerIterations");

        dict.readIfPresent("chargedSpecies", nCorrCharged_);
        dict.readIfPresent("neutralSpecies", nCorrNeutral_);
        dict.readIfPresent("electronTemperature", nCorrTe_);
        dict.readIfPresent("tolerance", innerTolerance_);
        dict.readIfPresent("reportUnconverged", reportUnconverged_);

        Info<< "Inner iterations: charged species " << nCorrCharged_
            << ", neutral species " << nCorrNeutral_
            << ", electron temperature " << nCorrTe_
            << ", tolerance " << innerTolerance_ << endl;
    }

    if (!found("slowSpeciesAcceleration"))
    {
        return;
    }

    const dictionary& dict = subDict("slowSpeciesAcceleration");

    accelerate_ = dict.lookupOrDefault<Switch>("active", true);

    if (!accelerate_)
    {
        return;
    }

    const wordList names(dict.lookup("species"));

    accSpecies_.setSize(names.size());

    forAll(names, k)
    {
        if (!species().contains(names[k]))
        {
            FatalIOErrorIn
            (
                "multiSpeciesPlasmaModel::readNumericalControls()",
                dict
            )   << "Unknown species " << names[k]
                << exit(FatalIOError);
        }

        accSpecies_[k] = species()[names[k]];

        if (accSpecies_[k] < activeSpecies_ || accSpecies_[k] == bIndex_)
        {
            FatalIOErrorIn
            (
                "multiSpeciesPlasmaModel::readNumericalControls()",
                dict
            )   << names[k] << " is a charged species or the background gas;"
                << " only neutral species transported by diffusion can be"
                << " accelerated" << exit(FatalIOError);
        }
    }

    accPeriod_ = readScalar(dict.lookup("period"));
    accFullCycles_ = readLabel(dict.lookup("fullCycles"));
    accAverageCycles_ = dict.lookupOrDefault<label>("averageCycles", 1);
    accDeltaT_ = readScalar(dict.lookup("deltaT"));
    accNSteps_ = readLabel(dict.lookup("nSteps"));
    accMaxChangeFactor_ = dict.lookupOrDefault<scalar>("maxChangeFactor", 2);
    accTolerance_ = dict.lookupOrDefault<scalar>("tolerance", 0);

    if
    (
        accPeriod_ <= 0
     || accDeltaT_ <= 0
     || accNSteps_ < 1
     || accAverageCycles_ < 1
     || accFullCycles_ < accAverageCycles_
     || accMaxChangeFactor_ <= 1
    )
    {
        FatalIOErrorIn
        (
            "multiSpeciesPlasmaModel::readNumericalControls()",
            dict
        )   << "Need period > 0, deltaT > 0, nSteps >= 1,"
            << " 1 <= averageCycles <= fullCycles and maxChangeFactor > 1"
            << exit(FatalIOError);
    }

    accSu_.setSize(names.size());
    accSp_.setSize(names.size());

    forAll(names, k)
    {
        accSu_.set
        (
            k,
            new volScalarField
            (
                IOobject
                (
                    "accSu_" + names[k],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        accSp_.set
        (
            k,
            new volScalarField
            (
                IOobject
                (
                    "accSp_" + names[k],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    Info<< "Slow species acceleration: " << names << " advanced by "
        << accNSteps_*accDeltaT_ << " s after every " << accFullCycles_
        << " periods of " << accPeriod_ << " s" << endl;
}


void Foam::multiSpeciesPlasmaModel::accelerateSlowSpecies
(
    psiChemistryModel& chemistry
)
{
    if (!accelerate_ || accConverged_)
    {
        return;
    }

    const scalar t = runTime_.value();
    const scalar dt = runTime_.deltaTValue();

    // Once per time step, whatever the number of outer correctors
    if (t <= accLastTime_)
    {
        return;
    }

    accLastTime_ = t;

    if (accBlockStart_ < -0.5*GREAT)
    {
        accBlockStart_ = t - dt;
    }

    const scalar blockLength = accFullCycles_*accPeriod_;
    const scalar timeInBlock = t - accBlockStart_;

    // Accumulate the sources over the last periods of the block. The
    // species equation is
    //     ddt(N) = laplacian(D, N) + Su + Sp*N
    // with Su = RR*A/W - dRRDi*N and Sp = dRRDi
    if (timeInBlock > blockLength - accAverageCycles_*accPeriod_ + 0.5*dt)
    {
        forAll(accSpecies_, k)
        {
            const label i = accSpecies_[k];

            const scalarField Sp(chemistry.dRRDi(i)().internalField());

            accSu_[k].internalField() +=
                dt
               *(
                    chemistry.RR(i)().internalField()
                   *plasmaConstants::A.value()/W(i)
                  - Sp*N_[i].internalField()
                );

            accSp_[k].internalField() += dt*Sp;
        }

        accAveragedTime_ += dt;
    }

    if (timeInBlock < blockLength - 0.5*dt)
    {
        return;
    }

    // Advance the slow species alone
    const dimensionedScalar rDeltaT
    (
        "rDeltaT",
        dimless/dimTime,
        1.0/accDeltaT_
    );

    scalar maxChange = 0;

    forAll(accSpecies_, k)
    {
        const label i = accSpecies_[k];

        volScalarField& Ni = N_[i];

        accSu_[k].internalField() /= accAveragedTime_;
        accSp_[k].internalField() /= accAveragedTime_;
        accSp_[k].correctBoundaryConditions();

        const scalarField N0(Ni.internalField());

        for (label stepI = 0; stepI < accNSteps_; stepI++)
        {
            const scalarField Nk(Ni.internalField());

            fvScalarMatrix NEqn
            (
                fvm::Sp(rDeltaT, Ni)
              - fvm::laplacian(D_[i], Ni, "laplacian(D,Nin)")
              - fvm::SuSp(accSp_[k], Ni)
            );

            NEqn.source() +=
                mesh_.V()*(rDeltaT.value()*Nk + accSu_[k].internalField());

            NEqn.solve(mesh_.solutionDict().solver("Nin"));

            Ni.max(1e4);
        }

        // Limit the change, since the plasma has not yet responded to it
        scalarField& NiI = Ni.internalField();

        NiI = min(max(NiI, N0/accMaxChangeFactor_), N0*accMaxChangeFactor_);

        Ni.correctBoundaryConditions();

        const scalar change =
            gMax(mag(NiI - N0))/max(gMax(N0), VSMALL);

        maxChange = max(maxChange, change);

        Info<< "Slow species acceleration: " << species()[i]
            << " advanced by " << accNSteps_*accDeltaT_
            << " s, maximum density " << gMax(N0) << " -> " << gMax(NiI)
            << ", largest change relative to the maximum " << change
            << endl;

        // The jump is not part of the next time derivative
        Ni.oldTime() == Ni;

        thermo_.composition().Y(i) =
            Ni*W(i)/thermo_.rho()/plasmaConstants::A;

        accSu_[k].internalField() = 0;
        accSp_[k].internalField() = 0;
    }

    updateChemistryCollFreq(chemistry);

    accAveragedTime_ = 0;
    accBlockStart_ = t;

    if (maxChange < accTolerance_)
    {
        accConverged_ = true;

        Info<< "Slow species acceleration: change below the tolerance "
            << accTolerance_ << ", no further advances" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiSpeciesPlasmaModel::reportInnerIterations
(
    const word& name,
    const label maxPasses,
    const scalar initialResidual
) const
{
    if (reportUnconverged_ && initialResidual >= innerTolerance_)
    {
        Info<< "Inner iterations of " << name << " not converged after "
            << maxPasses << " passes: initial residual of the last pass "
            << initialResidual << " (tolerance " << innerTolerance_ << ")"
            << endl;
    }
}


// ************************************************************************* //
