/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "zeroD.H"
#include "plasmaInterpolateXY.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::zeroD<ThermoType>::zeroD
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

    driveMode_("none"),

    driveTable_(),

    volume_(0),

    maxReducedField_(3000),

    powerField_(mesh_.nCells(), 0.0),

    wall_(false),

    radius_(0),

    length_(0),

    sigmaIn_(0),

    kRec_(0),

    wallIon_(),

    wallIonProduct_(),

    wallNeutral_(),

    wallGamma_(),

    wallDN_(),

    wallFromPool_(),

    wallProducts_(),

    wallPool_(-1),

    electronWallEnergyLoss_()
{
    if (found("zeroDWall"))
    {
        readWall(subDict("zeroDWall"));
    }

    if (found("zeroDDrive"))
    {
        const dictionary& driveDict = subDict("zeroDDrive");

        driveMode_ = word(driveDict.lookup("mode"));

        if
        (
            driveMode_ != "none"
         && driveMode_ != "reducedField"
         && driveMode_ != "electronTemperature"
         && driveMode_ != "power"
        )
        {
            FatalIOErrorIn("zeroD::zeroD", driveDict)
                << "unknown zeroDDrive mode " << driveMode_
                << ", use none, reducedField, electronTemperature or power"
                << exit(FatalIOError);
        }

        if (driveMode_ == "power")
        {
            if (driveDict.found("volume"))
            {
                volume_ = readScalar(driveDict.lookup("volume"));
            }
            else if (!wall_)
            {
                FatalIOErrorIn("zeroD::zeroD", driveDict)
                    << "the power drive needs a volume or zeroDWall"
                    << exit(FatalIOError);
            }

            maxReducedField_ =
                driveDict.lookupOrDefault<scalar>("maxReducedField", 3000);
        }

        if (driveDict.found("table"))
        {
            driveTable_ =
                List<Tuple2<scalar, scalar> >(driveDict.lookup("table"));
        }
        else if (driveMode_ != "none")
        {
            driveTable_.setSize(1);
            driveTable_[0] = Tuple2<scalar, scalar>
            (
                0,
                readScalar(driveDict.lookup("value"))
            );
        }

        if (driveMode_ != "none" && driveTable_.empty())
        {
            FatalIOErrorIn("zeroD::zeroD", driveDict)
                << "zeroDDrive needs a table or a value"
                << exit(FatalIOError);
        }

        Info<< "zeroD drive: " << driveMode_ << endl;
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::zeroD<ThermoType>::driveValue(const scalar t) const
{
    const List<Tuple2<scalar, scalar> >& tab = driveTable_;

    if (t <= tab[0].first())
    {
        return tab[0].second();
    }

    for (label k = 1; k < tab.size(); k++)
    {
        if (t <= tab[k].first())
        {
            const scalar w =
                (t - tab[k-1].first())
               /max(tab[k].first() - tab[k-1].first(), VSMALL);

            return (1 - w)*tab[k-1].second() + w*tab[k].second();
        }
    }

    return tab[tab.size()-1].second();
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::setDrivenTe()
{
    volScalarField& Te = thermo_.Te();

    Te == dimensionedScalar
    (
        "Te",
        Te.dimensions(),
        driveValue(mesh_.time().value())
    );

    updateTemperature();
}

template<class ThermoType>
void Foam::zeroD<ThermoType>::readWall(const dictionary& dict)
{
    const scalar pi = 3.14159265358979323846;

    wall_ = true;
    radius_ = readScalar(dict.lookup("radius"));
    length_ = readScalar(dict.lookup("length"));
    sigmaIn_ = readScalar(dict.lookup("ionNeutralCrossSection"));
    kRec_ = readScalar(dict.lookup("ionIonRecombination"));

    if (volume_ <= 0)
    {
        volume_ = pi*sqr(radius_)*length_;
    }

    const speciesTable& sp = thermo_.composition().species();

    List<Tuple2<word, word> > ions(dict.lookup("ions"));

    wallIon_.setSize(ions.size());
    wallIonProduct_.setSize(ions.size());

    forAll(ions, k)
    {
        wallIon_[k] = sp[ions[k].first()];
        wallIonProduct_[k] = sp[ions[k].second()];
    }

    const dictionary& nd = dict.subDict("neutrals");
    const wordList names(nd.toc());

    wallNeutral_.setSize(names.size());
    wallGamma_.setSize(names.size());
    wallDN_.setSize(names.size());
    wallFromPool_.setSize(names.size());
    wallProducts_.setSize(names.size());

    forAll(names, k)
    {
        const dictionary& d = nd.subDict(names[k]);

        wallNeutral_[k] = sp[names[k]];
        wallGamma_[k] = readScalar(d.lookup("gamma"));
        wallDN_[k] = readScalar(d.lookup("DN"));
        wallFromPool_[k] = d.lookupOrDefault<scalar>("fromPool", 0);

        if (d.lookupOrDefault<Switch>("pool", false))
        {
            if (wallPool_ >= 0)
            {
                FatalIOErrorIn("zeroD::readWall", d)
                    << "only one pool species is allowed"
                    << exit(FatalIOError);
            }
            wallPool_ = k;
        }

        List<Tuple2<word, scalar> > prod(d.lookup("products"));

        wallProducts_[k].setSize(prod.size());

        forAll(prod, j)
        {
            wallProducts_[k][j] =
                Tuple2<label, scalar>(sp[prod[j].first()], prod[j].second());
        }
    }

    forAll(wallNeutral_, k)
    {
        if (wallFromPool_[k] > 0 && wallPool_ < 0)
        {
            FatalIOErrorIn("zeroD::readWall", dict)
                << names[k] << " needs a pool species"
                << exit(FatalIOError);
        }
    }

    electronWallEnergyLoss_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "electronWallEnergyLoss",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    Info<< "zeroD walls: R " << radius_ << " m, L " << length_ << " m, "
        << wallIon_.size() << " ions, " << wallNeutral_.size()
        << " neutrals" << endl;
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::wallLosses(const scalar dt)
{
    const scalar pi = 3.14159265358979323846;
    const scalar kB = plasmaConstants::boltzC.value();
    const scalar Av = plasmaConstants::A.value();

    const scalar R = radius_;
    const scalar L = length_;
    const scalar V = pi*sqr(R)*L;
    const scalar Aw = 2*pi*R*(R + L);
    const scalar rLambda2 = sqr(pi/L) + sqr(2.405/R);

    const scalarField& Te = thermo_.Te().internalField();
    const scalarField& Tg = thermo_.T().internalField();
    const volScalarField rhoF(thermo_.rho());
    const scalarField& rho = rhoF.internalField();

    scalarField& eLoss = electronWallEnergyLoss_().internalField();

    const label nSp = species().size();

    forAll(Te, c)
    {
        scalar Ng = NG_.internalField()[c];

        if (Ng <= 0)
        {
            forAll(N_, i)
            {
                Ng += N_[i].internalField()[c];
            }
        }

        const scalar ne = N_[eIndex_].internalField()[c];

        // Positive and negative ions, and the mean mass of the positive ions
        scalar nPos = 0, nNeg = 0, nm = 0;

        for (label i = 0; i < nSp; i++)
        {
            if (i == eIndex_)
            {
                continue;
            }

            const scalar n = N_[i].internalField()[c];

            if (z(i) > 0)
            {
                nPos += n;
                nm += n*W(i)/Av;
            }
            else if (z(i) < 0)
            {
                nNeg += n;
            }
        }

        scalarField dN(nSp, 0.0);

        scalar eLost = 0;

        if (nPos > 0)
        {
            // Edge-to-centre ratios, ions at the gas temperature (eta = 1,
            // gamma_- = gamma_+ = Te/Tg)
            const scalar alpha = nNeg/max(ne, VSMALL);
            const scalar alpha0 = 1.5*alpha;
            const scalar lambdai = 1/(Ng*sigmaIn_);
            const scalar eta = 1;
            const scalar gam = Te[c]/Tg[c];
            const scalar mBar = nm/nPos;
            const scalar vi = sqrt(8*kB*Tg[c]/(pi*mBar));
            const scalar nStar = 15.0/56.0*sqr(eta)/(kRec_*lambdai)*vi;

            const scalar hc =
                nNeg > 0
              ? 1/(sqrt(gam) + sqrt(gam)*sqrt(nStar)*nPos/pow(nNeg, 1.5))
              : 0;

            const scalar hL = sqrt
            (
                sqr(0.86/sqrt(3 + eta*L/(2*lambdai))/(1 + alpha0)) + sqr(hc)
            );
            const scalar hR = sqrt
            (
                sqr(0.8/sqrt(4 + eta*R/lambdai)/(1 + alpha0)) + sqr(hc)
            );
            const scalar Aeff = 2*pi*(sqr(R)*hL + R*L*hR);

            forAll(wallIon_, k)
            {
                const label i = wallIon_[k];
                const scalar uB = sqrt(kB*Te[c]/(W(i)/Av));
                const scalar n = N_[i].internalField()[c];
                const scalar lost = n*(1 - exp(-uB*Aeff/V*dt));

                dN[i] -= lost;
                dN[wallIonProduct_[k]] += lost;
                eLost += lost;
            }

            eLost = min(eLost, ne);
            dN[eIndex_] -= eLost;
        }

        eLoss[c] = 2*kB*Te[c]*eLost/dt;

        // Neutrals
        scalarList lost(wallNeutral_.size(), 0.0);

        forAll(wallNeutral_, k)
        {
            const label i = wallNeutral_[k];
            const scalar g = wallGamma_[k];
            const scalar v = sqrt(8*kB*Tg[c]/(pi*W(i)/Av));
            const scalar D = wallDN_[k]/Ng;
            const scalar kw =
                1/(1/(rLambda2*D) + 2*V*(2 - g)/(Aw*v*g));

            lost[k] = N_[i].internalField()[c]*(1 - exp(-kw*dt));
        }

        // Share of the pool atoms
        scalar demand = 0;

        forAll(wallNeutral_, k)
        {
            demand += wallFromPool_[k]*lost[k];
        }

        const scalar poolAtoms = wallPool_ >= 0 ? lost[wallPool_] : 0;
        const scalar f = demand > 0 ? min(1.0, poolAtoms/demand) : 1.0;

        forAll(wallNeutral_, k)
        {
            scalar returned;

            if (k == wallPool_)
            {
                // The pool atoms left over return as the pool products
                returned = lost[k] - f*demand;
                dN[wallNeutral_[k]] -= lost[k];
            }
            else if (wallFromPool_[k] > 0)
            {
                returned = f*lost[k];
                dN[wallNeutral_[k]] -= returned;
            }
            else
            {
                returned = lost[k];
                dN[wallNeutral_[k]] -= returned;
            }

            forAll(wallProducts_[k], j)
            {
                dN[wallProducts_[k][j].first()] +=
                    returned*wallProducts_[k][j].second();
            }
        }

        // Densities and mass fractions; the background gas is set from the
        // others afterwards
        forAll(dN, i)
        {
            if (dN[i] != 0 && i != bIndex_)
            {
                volScalarField& yi = thermo_.composition().Y(i);

                N_[i].internalField()[c] =
                    max(N_[i].internalField()[c] + dN[i], 0.0);

                yi.internalField()[c] =
                    N_[i].internalField()[c]*W(i)/(Av*rho[c]);
            }
        }
    }

    forAll(N_, i)
    {
        N_[i].correctBoundaryConditions();
        thermo_.composition().Y(i).correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::zeroD<ThermoType>::input()
{
    // The species dictionaries (mobility models among them) are read here
    multiSpeciesPlasmaModel::input();

    if
    (
        (driveMode_ == "reducedField" || driveMode_ == "power")
     && mobilityModel_[eIndex_] != "eTemp"
     && mobilityModel_[eIndex_] != "constant"
    )
    {
        FatalErrorIn("zeroD::input")
            << "reducedField and power drives need the electron mobilityModel"
            << " eTemp or constant, not "
            << mobilityModel_[eIndex_]
            << exit(FatalError);
    }

    if (driveMode_ == "electronTemperature")
    {
        setDrivenTe();
    }
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::setDrivenField(volVectorField& E) const
{
    E == dimensionedVector("zero", E.dimensions(), vector::zero);

    if (driveMode_ == "reducedField")
    {
        // E/N in Td (1e-21 V m^2) times the gas density
        const scalar EN = 1e-21*driveValue(mesh_.time().value());

        E.internalField().replace(vector::X, EN*NG_.internalField());
    }
    else if (driveMode_ == "power")
    {
        // e Ne mu_e E^2 = P/V with Ne and mu_e(Te) of the start of the step
        const scalar PV = driveValue(mesh_.time().value())/volume_;
        const scalar e = plasmaConstants::eCharge.value();

        // Gas density: NG_ is set only by the first species update, so
        // before it the sum of the densities is taken
        scalarField Ng(NG_.internalField());

        if (gMin(Ng) <= 0)
        {
            Ng = 0.0;

            forAll(N_, i)
            {
                Ng += N_[i].internalField();
            }
        }

        scalarField mue(mesh_.nCells(), 0.0);

        if (mobilityModel_[eIndex_] == "eTemp")
        {
            mue = plasmaInterpolateXY
            (
                thermo_.Te().internalField(),
                graphMuData_[eIndex_].x(),
                graphMuData_[eIndex_].y()
            )/Ng;
        }
        else
        {
            mue = mu_[eIndex_].internalField();
        }

        const scalarField& Ne = N_[eIndex_].internalField();

        forAll(powerField_, c)
        {
            const scalar Emax = 1e-21*maxReducedField_*Ng[c];
            const scalar denom = e*Ne[c]*mag(mue[c]);

            powerField_[c] =
                denom > 0 ? min(sqrt(max(PV, 0.0)/denom), Emax) : Emax;

            if (PV <= 0)
            {
                powerField_[c] = 0;
            }
        }

        E.internalField().replace(vector::X, powerField_);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::zeroD<ThermoType>::correct
(
    psiChemistryModel& chemistry,
    const volVectorField& E,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{
    if (multiTimeStep)
    {
        FatalError << "0-D model not implemented for multiple time scales" << nl << abort(FatalError);
    }

    // The base class has evaluated the chemistry of the state at the start
    // of the step: Sy_ = RR, the rate averaged over the step with
    // plasmaChemistryModel temporal. Without transport the species change
    // by it alone, Y = Y0 + dt RR/rho
    const scalar dt = mesh_.time().deltaTValue();

    // Store the densities of the start of the step before they change: the
    // electron energy equation takes Ne.oldTime(), which would otherwise be
    // stored only at its first use, after this update
    forAll(N_, i)
    {
        N_[i].oldTime();
    }

    scalarField yt(mesh_.nCells(), 0.0);

    forAll(species(), i)
    {
        if (i != bIndex_ && speciesSolution_[i])
        {
            volScalarField& yi = thermo_.composition().Y(i);

            yi.internalField() =
                max
                (
                    yi.oldTime().internalField()
                  + dt*Sy_[i].internalField()/thermo_.rho()().internalField(),
                    0.0
                );

            yi.correctBoundaryConditions();

            yt += yi.internalField();

            N_[i] == thermo_.rho()*yi*plasmaConstants::A/W(i);
        }
        else if (i != bIndex_)
        {
            yt += thermo_.composition().Y(i).internalField();
        }
    }

    // Losses to the walls (zeroDWall), with the densities of the start of
    // the step for the electrons and Te
    if (wall_)
    {
        wallLosses(dt);

        yt = 0.0;

        forAll(species(), i)
        {
            if (i != bIndex_)
            {
                yt += thermo_.composition().Y(i).internalField();
            }
        }
    }

    volScalarField& yBgas = thermo_.composition().Y()[bIndex_];

    yBgas.internalField() = 1 - yt;

    yBgas.correctBoundaryConditions();

    N_[bIndex_] == thermo_.rho()*yBgas*plasmaConstants::A/W(bIndex_);

    // Electron mobility at the current Te, also needed by the electron
    // conductivity and the collision frequency
    volScalarField& mue = mu_[eIndex_];

    if (mobilityModel_[eIndex_] == "eTemp")
    {
        mue.field() = (1/NG_.internalField())*plasmaInterpolateXY
        (
            thermo_.Te().internalField(),
            graphMuData_[eIndex_].x(),
            graphMuData_[eIndex_].y()
        );
    }

    // Electron flux in the driving field, whose work J.E heats the
    // electrons in the electron energy equation
    if (driveMode_ == "reducedField")
    {
        J_[eIndex_] == N_[eIndex_]*mue*z_[eIndex_]*E;
    }
    else if (driveMode_ == "power")
    {
        // Work J.E of the step equal to the power density P/V, or the
        // flux of the electrons in the limited field when they cannot take
        // that power
        const scalar PV = driveValue(mesh_.time().value())/volume_;
        const scalar e = plasmaConstants::eCharge.value();

        J_[eIndex_] == N_[eIndex_]*mue*z_[eIndex_]*E;

        vectorField& Je = J_[eIndex_].internalField();

        forAll(Je, c)
        {
            const scalar Ex = E.internalField()[c].x();

            if (Ex > 0)
            {
                Je[c].x() = -min(PV/(e*Ex), mag(Je[c].x()));
            }
        }
    }
    else
    {
        J_[eIndex_] == dimensionedVector
        (
            "zero",
            J_[eIndex_].dimensions(),
            vector::zero
        );
    }

    // Prescribed Te of the next step: the chemistry evaluated at the end of
    // this step uses it
    if (driveMode_ == "electronTemperature")
    {
        setDrivenTe();
    }

    return 0;
}

template<class ThermoType>
bool Foam::zeroD<ThermoType>::read()
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
