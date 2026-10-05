/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Description
    Numerical controls of multiSpeciesPlasmaModel read from plasmaProperties:

    Limits. densityFloor, TeMin and TeMax are the defaults for the number
    densities of the transported species and for the electron temperature
    (values shown; no upper limit on the electron temperature unless TeMax
    is given). A sub-dictionary named after a field sets a floor (min)
    and/or a ceiling (max) for that field alone and takes precedence:
    N_<specie> for the number densities, Te, Tion and T. The dictionary is
    read again when plasmaProperties is modified during the run.
    \verbatim
    limits
    {
        densityFloor    1e4;    // [1/m3]
        TeMin           300;    // [K]
        TeMax           1e6;    // [K], default: none

        N_electron      { min 1e10; }
        N_Arm           { min 1e12; max 1e20; }
        Te              { min 300; max 1.2e5; }
    }
    \endverbatim

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
        startTime       0;      // no advances before this time [s], e.g.
                                // until the discharge has ignited
        fullCycles      5;      // periods of full simulation per block
        averageCycles   1;      // periods at the end of a block over which
                                // the sources are averaged
        deltaT          1e-6;   // time step of the advance [s]
        nSteps          100;    // steps per advance
        steadyState     no;     // yes: solve the steady state of the
                                // averaged equation instead; deltaT and
                                // nSteps are then optional and used only
                                // where there is no steady state
        relaxation      1;      // fraction of the change that is applied
        maxChangeFactor 2;      // limit on the density change per advance
        tolerance       1e-3;   // advances stop below this relative change
    }
    \endverbatim

    After each block the listed species alone are advanced by
    nSteps*deltaT with the period-averaged chemistry source, linearised in
    the species' own density, and diffusion; the charged species, the
    electron temperature and the potential are left as they are and adjust
    during the next block.

    With steadyState the averaged equation
        - laplacian(D, N) - Sp*N = Su
    is solved once per block instead. Sources that are nonlinear in the
    species' own density (e.g. metastable pooling) are linearised about the
    density of the block, so the blocks act as the nonlinear iteration; use
    relaxation and maxChangeFactor to keep the changes moderate. A steady
    state requires a net loss (Sp < 0) in every cell; otherwise the time
    step advance is used if deltaT and nSteps are given, and the advance is
    skipped if not.

\*---------------------------------------------------------------------------*/

#include "multiSpeciesPlasmaModel.H"
#include "snGradScheme.H"
#include "zeroGradientFvPatchFields.H"
#include "fvm.H"
#include "OSspecific.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::multiSpeciesPlasmaModel::driftDiffusionTerms
(
    const label i,
    const volVectorField& F,
    const volScalarField& D,
    const volScalarField& Ni
) const
{
    const surfaceScalarField phi(fvc::interpolate(F) & mesh_.Sf());

    if (scharfetterGummel_[i])
    {
        return Foam::scharfetterGummel(phi, D, Ni);
    }

    return
        fvm::div(phi, Ni, "div(F,Ni)")
      - fvm::laplacian(D, Ni, "laplacian(D,Ni)");
}


void Foam::multiSpeciesPlasmaModel::correctDriftDiffusionFlux
(
    const label i,
    const volVectorField& F,
    const volScalarField& D,
    const volScalarField& Ni
)
{
    const surfaceScalarField phi(fvc::interpolate(F) & mesh_.Sf());

    // The flux as discretised in the equation of the species
    tmp<surfaceScalarField> tflux;

    if (scharfetterGummel_[i])
    {
        tflux = Foam::scharfetterGummelFlux(phi, D, Ni);
    }
    else
    {
        // Interpolation and surface-normal gradient schemes of the
        // Laplacian term, "Gauss <interpolation> <snGrad>"
        ITstream& is = mesh_.schemesDict().laplacianScheme("laplacian(D,Ni)");

        const word gauss(is);

        tmp<surfaceInterpolationScheme<scalar> > tinterpolation
        (
            surfaceInterpolationScheme<scalar>::New(mesh_, is)
        );

        tmp<fv::snGradScheme<scalar> > tsnGrad
        (
            fv::snGradScheme<scalar>::New(mesh_, is)
        );

        tflux =
            fvc::flux(phi, Ni, "div(F,Ni)")
          - tinterpolation().interpolate(D)*mesh_.magSf()
           *tsnGrad().snGrad(Ni);
    }

    // Not registered: it is set again after every solution of the species,
    // also after the mesh has changed
    faceFlux_.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "faceFlux_" + species()[i],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            tflux
        )
    );

    // The boundary values of J, which give the currents to the walls, are
    // kept
    J_[i].internalField() = fvc::reconstruct(faceFlux_[i])().internalField();
}


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

    // Flux scheme: for all species, and optionally per species
    {
        const word allScheme(lookupOrDefault<word>("fluxScheme", "fvSchemes"));

        Info<< "Drift-diffusion fluxes: " << allScheme;

        forAll(species(), i)
        {
            word scheme(allScheme);

            if (isDict(species()[i]))
            {
                subDict(species()[i]).readIfPresent("fluxScheme", scheme);
            }

            if (scheme != "fvSchemes" && scheme != "scharfetterGummel")
            {
                FatalIOErrorIn
                (
                    "multiSpeciesPlasmaModel::readNumericalControls()",
                    *this
                )   << "Unknown fluxScheme " << scheme
                    << "; valid entries are fvSchemes and scharfetterGummel"
                    << exit(FatalIOError);
            }

            scharfetterGummel_[i] = (scheme == "scharfetterGummel");

            if (i < activeSpecies_ && scheme != allScheme)
            {
                Info<< ", " << species()[i] << " " << scheme;
            }
        }

        Info<< endl;
    }

    limitsTime_ = lastModified(filePath());

    readLimits(*this);

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
    accStartTime_ = dict.lookupOrDefault<scalar>("startTime", 0);
    accFullCycles_ = readLabel(dict.lookup("fullCycles"));
    accAverageCycles_ = dict.lookupOrDefault<label>("averageCycles", 1);
    accSteady_ = dict.lookupOrDefault<Switch>("steadyState", false);
    accRelaxation_ = dict.lookupOrDefault<scalar>("relaxation", 1);

    if (accSteady_)
    {
        // Only used where there is no steady state
        accDeltaT_ = dict.lookupOrDefault<scalar>("deltaT", 0);
        accNSteps_ = dict.lookupOrDefault<label>("nSteps", 0);
    }
    else
    {
        accDeltaT_ = readScalar(dict.lookup("deltaT"));
        accNSteps_ = readLabel(dict.lookup("nSteps"));
    }
    accMaxChangeFactor_ = dict.lookupOrDefault<scalar>("maxChangeFactor", 2);
    accTolerance_ = dict.lookupOrDefault<scalar>("tolerance", 0);

    if
    (
        accPeriod_ <= 0
     || (!accSteady_ && (accDeltaT_ <= 0 || accNSteps_ < 1))
     || accRelaxation_ <= 0
     || accRelaxation_ > 1
     || accAverageCycles_ < 1
     || accFullCycles_ < accAverageCycles_
     || accMaxChangeFactor_ <= 1
    )
    {
        FatalIOErrorIn
        (
            "multiSpeciesPlasmaModel::readNumericalControls()",
            dict
        )   << "Need period > 0, deltaT > 0 and nSteps >= 1 (unless"
            << " steadyState), 0 < relaxation <= 1,"
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

    Info<< "Slow species acceleration:";

    forAll(names, k)
    {
        Info<< " " << names[k];
    }

    if (accSteady_)
    {
        Info<< " solved for their steady state";
    }
    else
    {
        Info<< " advanced by " << accNSteps_*accDeltaT_ << " s";
    }

    Info<< " after every " << accFullCycles_ << " periods of " << accPeriod_
        << " s";

    if (accStartTime_ > 0)
    {
        Info<< ", from time " << accStartTime_ << " s";
    }

    Info<< endl;
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

    // The first block starts at startTime
    if (t < accStartTime_ + 0.5*dt)
    {
        return;
    }

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
        accDeltaT_ > 0 ? 1.0/accDeltaT_ : 0.0
    );

    scalar maxChange = 0;

    forAll(accSpecies_, k)
    {
        const label i = accSpecies_[k];

        volScalarField& Ni = N_[i];

        accSu_[k].internalField() /= accAveragedTime_;
        accSp_[k].internalField() /= accAveragedTime_;
        accSp_[k].correctBoundaryConditions();

        // The source is Su + Sp*N. fvm::SuSp(-Sp, N) on the left-hand side
        // makes a loss (Sp < 0) implicit and a gain explicit

        const scalarField N0(Ni.internalField());

        // A steady state needs a net loss everywhere
        bool steady = accSteady_;

        if (steady && gMax(accSp_[k].internalField()) >= 0)
        {
            steady = false;

            WarningIn("multiSpeciesPlasmaModel::accelerateSlowSpecies(...)")
                << species()[i] << " has no net loss in part of the domain"
                << " (largest averaged dRRDi "
                << gMax(accSp_[k].internalField()) << "): no steady state, "
                << (accNSteps_ > 0 ? "advancing in time" : "advance skipped")
                << endl;
        }

        if (steady)
        {
            fvScalarMatrix NEqn
            (
              - fvm::laplacian(D_[i], Ni, "laplacian(D,Nin)")
              + fvm::SuSp(-accSp_[k], Ni)
            );

            NEqn.source() += mesh_.V()*accSu_[k].internalField();

            NEqn.solve(mesh_.solutionDict().solver("Nin"));

            limitField(Ni, densityFloor_);
        }

        for
        (
            label stepI = 0;
            !steady && accDeltaT_ > 0 && stepI < accNSteps_;
            stepI++
        )
        {
            const scalarField Nk(Ni.internalField());

            fvScalarMatrix NEqn
            (
                fvm::Sp(rDeltaT, Ni)
              - fvm::laplacian(D_[i], Ni, "laplacian(D,Nin)")
              + fvm::SuSp(-accSp_[k], Ni)
            );

            NEqn.source() +=
                mesh_.V()*(rDeltaT.value()*Nk + accSu_[k].internalField());

            NEqn.solve(mesh_.solutionDict().solver("Nin"));

            limitField(Ni, densityFloor_);
        }

        // Limit the change, since the plasma has not yet responded to it
        scalarField& NiI = Ni.internalField();

        NiI = N0 + accRelaxation_*(NiI - N0);

        NiI = min(max(NiI, N0/accMaxChangeFactor_), N0*accMaxChangeFactor_);

        Ni.correctBoundaryConditions();

        const scalar change =
            gMax(mag(NiI - N0))/max(gMax(N0), VSMALL);

        maxChange = max(maxChange, change);

        Info<< "Slow species acceleration: " << species()[i];

        if (steady)
        {
            Info<< " solved for steady state";
        }
        else
        {
            Info<< " advanced by " << accNSteps_*accDeltaT_ << " s";
        }

        Info<< ", maximum density " << gMax(N0) << " -> " << gMax(NiI)
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


void Foam::multiSpeciesPlasmaModel::readLimitsIfModified()
{
    if
    (
       !runTime_.controlDict().lookupOrDefault<Switch>
        (
            "runTimeModifiable",
            true
        )
    )
    {
        return;
    }

    const time_t fileTime = lastModified(filePath());

    if (fileTime <= limitsTime_)
    {
        return;
    }

    limitsTime_ = fileTime;

    // plasmaProperties may be registered more than once, so that this
    // object is not told about the modification: read the file directly
    IOdictionary properties
    (
        IOobject
        (
            name(),
            instance(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Info<< "Reading the limits again from the modified " << name() << endl;

    readLimits(properties);
}


void Foam::multiSpeciesPlasmaModel::readLimits(const dictionary& properties)
{
    fieldMin_.clear();
    fieldMax_.clear();

    if (!properties.found("limits"))
    {
        return;
    }

    const dictionary& dict = properties.subDict("limits");

    dict.readIfPresent("densityFloor", densityFloor_);
    dict.readIfPresent("TeMin", TeMin_);
    dict.readIfPresent("TeMax", TeMax_);

    if (densityFloor_ < 0 || TeMin_ <= 0 || TeMax_ <= TeMin_)
    {
        FatalIOErrorIn("multiSpeciesPlasmaModel::readLimits()", dict)
            << "Need densityFloor >= 0 and 0 < TeMin < TeMax"
            << exit(FatalIOError);
    }

    Info<< "Limits: densities not below " << densityFloor_
        << " 1/m3, electron temperature between " << TeMin_ << " and ";

    if (TeMax_ < 0.5*GREAT)
    {
        Info<< TeMax_ << " K" << endl;
    }
    else
    {
        Info<< "unlimited" << endl;
    }

    // Limits of individual fields: a sub-dictionary named after the field
    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            continue;
        }

        const word& fieldName = iter().keyword();
        const dictionary& fieldDict = iter().dict();

        if (fieldDict.found("min"))
        {
            fieldMin_.insert(fieldName, readScalar(fieldDict.lookup("min")));
        }

        if (fieldDict.found("max"))
        {
            fieldMax_.insert(fieldName, readScalar(fieldDict.lookup("max")));
        }

        if
        (
            fieldMin_.found(fieldName)
         && fieldMax_.found(fieldName)
         && fieldMax_[fieldName] <= fieldMin_[fieldName]
        )
        {
            FatalIOErrorIn("multiSpeciesPlasmaModel::readLimits()", fieldDict)
                << "max must be above min for " << fieldName
                << exit(FatalIOError);
        }

        Info<< "Limits: " << fieldName;

        if (fieldMin_.found(fieldName))
        {
            Info<< " not below " << fieldMin_[fieldName];
        }

        if (fieldMax_.found(fieldName))
        {
            Info<< " not above " << fieldMax_[fieldName];
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiSpeciesPlasmaModel::limitField
(
    volScalarField& field,
    const scalar defaultMin,
    const scalar defaultMax
) const
{
    const word& name = field.name();

    const scalar lower =
        fieldMin_.found(name) ? fieldMin_[name] : defaultMin;

    const scalar upper =
        fieldMax_.found(name) ? fieldMax_[name] : defaultMax;

    if (lower > -0.5*GREAT)
    {
        field.max(lower);
    }

    if (upper < 0.5*GREAT)
    {
        field.min(upper);
    }
}



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
