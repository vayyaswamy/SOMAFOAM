/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "zeroD.H"
#include "plasmaInterpolateXY.H"
#include "simpleMatrix.H"

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

    electronWallEnergyLoss_(),

    coupled_(false),

    debugCharge_(false),

    relTol_(1e-5),

    absTol_(1e6),

    hEst_(-1),

    chem_(NULL),

    cell_(0),

    eKind_()
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

        coupled_ = driveDict.lookupOrDefault<Switch>("coupled", false);
        relTol_ = driveDict.lookupOrDefault<scalar>("relTol", 1e-5);
        debugCharge_ =
            driveDict.lookupOrDefault<Switch>("debugCharge", false);
        absTol_ = driveDict.lookupOrDefault<scalar>("absTol", 1e6);

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

        Info<< "zeroD drive: " << driveMode_
            << (coupled_ ? ", coupled implicit step" : "")
            << endl;
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
void Foam::zeroD<ThermoType>::wallFrequencies
(
    const scalarField& n,
    const scalar Te,
    const scalar Tg,
    const scalar Ng,
    scalarList& kIon,
    scalarList& kNeutral
) const
{
    const scalar pi = 3.14159265358979323846;
    const scalar kB = plasmaConstants::boltzC.value();
    const scalar Av = plasmaConstants::A.value();

    const scalar R = radius_;
    const scalar L = length_;
    const scalar V = pi*sqr(R)*L;
    const scalar Aw = 2*pi*R*(R + L);
    const scalar rLambda2 = sqr(pi/L) + sqr(2.405/R);

    const scalar ne = n[eIndex_];

    // Positive and negative ions, and the mean mass of the positive ions
    scalar nPos = 0, nNeg = 0, nm = 0;

    forAll(n, i)
    {
        if (i == eIndex_)
        {
            continue;
        }

        if (z(i) > 0)
        {
            nPos += n[i];
            nm += n[i]*W(i)/Av;
        }
        else if (z(i) < 0)
        {
            nNeg += n[i];
        }
    }

    kIon.setSize(wallIon_.size());
    kIon = 0.0;

    if (nPos > 0)
    {
        // Edge-to-centre ratios, ions at the gas temperature (eta = 1,
        // gamma_- = gamma_+ = Te/Tg)
        const scalar alpha = nNeg/max(ne, VSMALL);
        const scalar alpha0 = 1.5*alpha;
        const scalar lambdai = 1/(Ng*sigmaIn_);
        const scalar eta = 1;
        const scalar gam = Te/Tg;
        const scalar mBar = nm/nPos;
        const scalar vi = sqrt(8*kB*Tg/(pi*mBar));
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
            kIon[k] = sqrt(kB*Te/(W(wallIon_[k])/Av))*Aeff/V;
        }
    }

    kNeutral.setSize(wallNeutral_.size());

    forAll(wallNeutral_, k)
    {
        const scalar g = wallGamma_[k];
        const scalar v = sqrt(8*kB*Tg/(pi*W(wallNeutral_[k])/Av));
        const scalar D = wallDN_[k]/Ng;

        kNeutral[k] = 1/(1/(rLambda2*D) + 2*V*(2 - g)/(Aw*v*g));
    }
}


template<class ThermoType>
Foam::scalar Foam::zeroD<ThermoType>::wallTransfer
(
    const scalarList& lostIon,
    const scalarList& lostNeutral,
    const scalar ne,
    scalarField& dN
) const
{
    // Ions return as their neutrals, with one electron each
    scalar eLost = 0;

    forAll(wallIon_, k)
    {
        dN[wallIon_[k]] -= lostIon[k];
        dN[wallIonProduct_[k]] += lostIon[k];
        eLost += lostIon[k];
    }

    eLost = min(eLost, ne);
    dN[eIndex_] -= eLost;

    // Share of the pool atoms
    scalar demand = 0;

    forAll(wallNeutral_, k)
    {
        demand += wallFromPool_[k]*lostNeutral[k];
    }

    const scalar poolAtoms = wallPool_ >= 0 ? lostNeutral[wallPool_] : 0;
    const scalar f = demand > 0 ? min(1.0, poolAtoms/demand) : 1.0;

    forAll(wallNeutral_, k)
    {
        scalar returned;

        if (k == wallPool_)
        {
            // The pool atoms left over return as the pool products
            returned = lostNeutral[k] - f*demand;
            dN[wallNeutral_[k]] -= lostNeutral[k];
        }
        else if (wallFromPool_[k] > 0)
        {
            returned = f*lostNeutral[k];
            dN[wallNeutral_[k]] -= returned;
        }
        else
        {
            returned = lostNeutral[k];
            dN[wallNeutral_[k]] -= returned;
        }

        forAll(wallProducts_[k], j)
        {
            dN[wallProducts_[k][j].first()] +=
                returned*wallProducts_[k][j].second();
        }
    }

    return eLost;
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::wallLosses(const scalar dt)
{
    const scalar kB = plasmaConstants::boltzC.value();
    const scalar Av = plasmaConstants::A.value();

    const scalarField& Te = thermo_.Te().internalField();
    const scalarField& Tg = thermo_.T().internalField();
    const volScalarField rhoF(thermo_.rho());
    const scalarField& rho = rhoF.internalField();

    scalarField& eLoss = electronWallEnergyLoss_().internalField();

    const label nSp = species().size();

    forAll(Te, c)
    {
        scalarField n(nSp);
        scalar Ng = NG_.internalField()[c];

        forAll(n, i)
        {
            n[i] = N_[i].internalField()[c];
        }

        if (Ng <= 0)
        {
            Ng = sum(n);
        }

        scalarList kIon, kNeutral;
        wallFrequencies(n, Te[c], Tg[c], Ng, kIon, kNeutral);

        // Amounts lost over the step, exact for the densities of its start
        scalarList lostIon(wallIon_.size()), lostNeutral(wallNeutral_.size());

        forAll(wallIon_, k)
        {
            lostIon[k] = n[wallIon_[k]]*(1 - exp(-kIon[k]*dt));
        }

        forAll(wallNeutral_, k)
        {
            lostNeutral[k] = n[wallNeutral_[k]]*(1 - exp(-kNeutral[k]*dt));
        }

        scalarField dN(nSp, 0.0);

        const scalar eLost = wallTransfer(lostIon, lostNeutral, n[eIndex_], dN);

        eLoss[c] = 2*kB*Te[c]*eLost/dt;

        // Densities and mass fractions; the background gas is set from the
        // others afterwards
        forAll(dN, i)
        {
            if (dN[i] != 0 && i != bIndex_)
            {
                volScalarField& yi = thermo_.composition().Y(i);

                N_[i].internalField()[c] = max(n[i] + dN[i], 0.0);

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


template<class ThermoType>
Foam::scalar Foam::zeroD<ThermoType>::coupledTe
(
    const scalar t,
    const scalarField& y
) const
{
    const label nSp = species().size();

    scalar Te;

    if (driveMode_ == "electronTemperature")
    {
        Te = driveValue(t);
    }
    else
    {
        // y[nSp] = 3/2 Ne Te in eV/m3
        const scalar EV_K =
            plasmaConstants::eCharge.value()/plasmaConstants::boltzC.value();

        Te = max(y[nSp], 0.0)/(1.5*max(y[eIndex_], VSMALL))*EV_K;
    }

    return min(max(Te, TeMin_), TeMax_);
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::coupledDerivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    const ODEChemistryModel<psiChemistryModel, ThermoType>& chem = *chem_;
    const PtrList<Reaction<ThermoType> >& reactions = chem.reactions();

    const scalar Av = plasmaConstants::A.value();
    const scalar e = plasmaConstants::eCharge.value();
    // Boltzmann constant of the elastic loss in ODEChemistryModel
    const scalar kBel = 1.38064852E-23;
    const scalar kB = plasmaConstants::boltzC.value();

    const label c = cell_;
    const label nSp = species().size();

    const scalar Tg = thermo_.T().internalField()[c];
    const scalar p = thermo_.p().internalField()[c];
    const scalar Te = coupledTe(t, y);

    // Densities [1/m3] and concentrations [kmol/m3]
    scalarField n(nSp);
    scalarField conc(nSp);

    forAll(n, i)
    {
        n[i] = max(y[i], 0.0);
        conc[i] = n[i]/Av;
    }

    scalar Ng = NG_.internalField()[c];

    if (Ng <= 0)
    {
        Ng = sum(n);
    }

    dydt = 0.0;

    if (eKind_.size() != reactions.size())
    {
        eKind_.setSize(reactions.size());

        forAll(reactions, r)
        {
            const word& type = reactions[r].type();

            eKind_[r] =
                (
                    type == "electronImpactInelasticArrheniusReaction"
                 || type == "electronImpactInelasticTabularReaction"
                )
              ? 1
              : (
                    type == "electronImpactElasticArrheniusReaction"
                 || type == "electronImpactElasticTabularReaction"
                )
              ? 2
              : 0;
        }
    }

    // Reactions and their electron energy loss [W/m3]
    const scalar We = W(eIndex_);
    scalar loss = 0;

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    forAll(reactions, r)
    {
        const Reaction<ThermoType>& R = reactions[r];

        const scalar om = chem.omega
        (
            R, conc, Tg, p, Te, Tg, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            dydt[R.lhs()[s].index] -= R.lhs()[s].stoichCoeff*om*Av;
        }

        forAll(R.rhs(), s)
        {
            dydt[R.rhs()[s].index] += R.rhs()[s].stoichCoeff*om*Av;
        }

        if (eKind_[r] == 1)
        {
            loss += R.lhs()[0].stoichCoeff*om*Av*R.deltaE();
        }
        else if (eKind_[r] == 2)
        {
            const scalar M = W(R.lhs()[1].index);

            loss += 3*om*Av*(We/M)*kBel*(Te - Tg);
        }
    }

    // Walls
    scalar eLost = 0;

    if (wall_)
    {
        scalarList kIon, kNeutral;
        wallFrequencies(n, Te, Tg, Ng, kIon, kNeutral);

        scalarList rIon(wallIon_.size()), rNeutral(wallNeutral_.size());

        forAll(wallIon_, k)
        {
            rIon[k] = kIon[k]*n[wallIon_[k]];
        }

        forAll(wallNeutral_, k)
        {
            rNeutral[k] = kNeutral[k]*n[wallNeutral_[k]];
        }

        scalarField dN(nSp, 0.0);

        // Rates: the electron loss is not limited by ne here
        eLost = wallTransfer(rIon, rNeutral, VGREAT, dN);

        forAll(dN, i)
        {
            dydt[i] += dN[i];
        }
    }

    // Electron energy density [eV/m3]
    if (driveMode_ != "electronTemperature")
    {
        const scalar Ne = n[eIndex_];

        scalar mue = 0;

        if (mobilityModel_[eIndex_] == "eTemp")
        {
            const scalarField muT
            (
                plasmaInterpolateXY
                (
                    scalarField(1, Te),
                    graphMuData_[eIndex_].x(),
                    graphMuData_[eIndex_].y()
                )
            );

            mue = muT[0]/Ng;
        }
        else
        {
            mue = mu_[eIndex_].internalField()[c];
        }

        scalar heat = 0;

        if (driveMode_ == "power")
        {
            const scalar Emax = 1e-21*maxReducedField_*Ng;

            heat = min
            (
                max(driveValue(t), 0.0)/volume_,
                e*mag(mue)*Ne*sqr(Emax)
            );
        }
        else if (driveMode_ == "reducedField")
        {
            heat = e*mag(mue)*Ne*sqr(1e-21*driveValue(t)*Ng);
        }

        dydt[nSp] = (heat - loss - 2*kB*Te*eLost)/e;
    }
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::coupledJacobian
(
    const scalar t,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    const label n = y.size();

    scalarField f0(n), f1(n);
    coupledDerivatives(t, y, f0);

    scalarField yp(y);

    for (label j = 0; j < n; j++)
    {
        const scalar h = 1e-7*max(mag(y[j]), absTol_);

        yp[j] = y[j] + h;
        coupledDerivatives(t, yp, f1);
        yp[j] = y[j];

        for (label i = 0; i < n; i++)
        {
            dfdy[i][j] = (f1[i] - f0[i])/h;
        }
    }

    // Time dependence of the drive
    const scalar ht = max(1e-6*max(mag(t), VSMALL), 1e-15);
    coupledDerivatives(t + ht, y, f1);
    dfdt = (f1 - f0)/ht;
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::rosenbrockStep
(
    scalar& x,
    scalarField& y,
    const scalarField& dydx,
    const scalarField& yScale,
    const scalar hTry,
    scalar& hDid,
    scalar& hNext
) const
{
    // Kaps-Rentrop method with Shampine's coefficients (Numerical Recipes,
    // stiff). The third and fourth stages take the derivative evaluated at
    // the third stage.
    static const scalar
        safety = 0.9, grow = 1.5, pgrow = -0.25,
        shrink = 0.5, pshrink = -1.0/3.0, errcon = 0.1296,
        gamma = 1.0/2.0,
        a21 = 2.0, a31 = 48.0/25.0, a32 = 6.0/25.0,
        c21 = -8.0, c31 = 372.0/25.0, c32 = 12.0/5.0,
        c41 = -112.0/125.0, c42 = -54.0/125.0, c43 = -2.0/5.0,
        b1 = 19.0/9.0, b2 = 1.0/2.0, b3 = 25.0/108.0, b4 = 125.0/108.0,
        e1 = 17.0/54.0, e2 = 7.0/36.0, e3 = 0.0, e4 = 125.0/108.0,
        c1X = 1.0/2.0, c2X = -3.0/2.0, c3X = 121.0/50.0, c4X = 29.0/250.0,
        a2X = 1.0, a3X = 3.0/5.0;

    const label n = y.size();
    const scalar x0 = x;
    const scalarField y0(y);

    scalarField dfdx(n);
    scalarSquareMatrix dfdy(n, 0.0);
    coupledJacobian(x0, y0, dfdx, dfdy);

    scalarSquareMatrix a(n, 0.0);
    labelList pivot(n);
    scalarField g1(n), g2(n), g3(n), g4(n), f(n), yErr(n);

    scalar h = hTry;

    for (label jtry = 0; jtry < 40; jtry++)
    {
        for (label i = 0; i < n; i++)
        {
            for (label j = 0; j < n; j++)
            {
                a[i][j] = -dfdy[i][j];
            }

            a[i][i] += 1.0/(gamma*h);
        }

        simpleMatrix<scalar>::LUDecompose(a, pivot);

        for (label i = 0; i < n; i++)
        {
            g1[i] = dydx[i] + h*c1X*dfdx[i];
        }

        simpleMatrix<scalar>::LUBacksubstitute(a, pivot, g1);

        for (label i = 0; i < n; i++)
        {
            y[i] = y0[i] + a21*g1[i];
        }

        coupledDerivatives(x0 + a2X*h, y, f);

        for (label i = 0; i < n; i++)
        {
            g2[i] = f[i] + h*c2X*dfdx[i] + c21*g1[i]/h;
        }

        simpleMatrix<scalar>::LUBacksubstitute(a, pivot, g2);

        for (label i = 0; i < n; i++)
        {
            y[i] = y0[i] + a31*g1[i] + a32*g2[i];
        }

        coupledDerivatives(x0 + a3X*h, y, f);

        for (label i = 0; i < n; i++)
        {
            g3[i] = f[i] + h*c3X*dfdx[i] + (c31*g1[i] + c32*g2[i])/h;
        }

        simpleMatrix<scalar>::LUBacksubstitute(a, pivot, g3);

        for (label i = 0; i < n; i++)
        {
            g4[i] =
                f[i] + h*c4X*dfdx[i]
              + (c41*g1[i] + c42*g2[i] + c43*g3[i])/h;
        }

        simpleMatrix<scalar>::LUBacksubstitute(a, pivot, g4);

        scalar maxErr = 0;

        for (label i = 0; i < n; i++)
        {
            y[i] = y0[i] + b1*g1[i] + b2*g2[i] + b3*g3[i] + b4*g4[i];
            yErr[i] = e1*g1[i] + e2*g2[i] + e3*g3[i] + e4*g4[i];
            maxErr = max(maxErr, mag(yErr[i]/yScale[i]));
        }

        maxErr /= relTol_;

        if (maxErr <= 1)
        {
            x = x0 + h;
            hDid = h;
            hNext = maxErr > errcon ? safety*h*pow(maxErr, pgrow) : grow*h;
            return;
        }

        hNext = safety*h*pow(maxErr, pshrink);
        h = max(hNext, shrink*h);
    }

    FatalErrorIn("zeroD::rosenbrockStep")
        << "no acceptable step after 40 tries at t = " << x0
        << exit(FatalError);
}


template<class ThermoType>
void Foam::zeroD<ThermoType>::coupledStep
(
    psiChemistryModel& chemistry,
    const scalar dt
)
{
    chem_ =
        dynamic_cast<const ODEChemistryModel<psiChemistryModel, ThermoType>*>
        (
            &chemistry
        );

    if (!chem_)
    {
        FatalErrorIn("zeroD::coupledStep")
            << "the coupled step needs an ODEChemistryModel"
            << exit(FatalError);
    }

    const scalar kB = plasmaConstants::boltzC.value();
    const scalar e = plasmaConstants::eCharge.value();
    const scalar Av = plasmaConstants::A.value();
    const label nSp = species().size();

    // The step ends at the current time
    const scalar t1 = mesh_.time().value();
    const scalar t0 = t1 - dt;

    volScalarField& TeF = thermo_.Te();
    const volScalarField rhoF(thermo_.rho());
    const scalarField& rho = rhoF.internalField();

    scalarField y(nSp + 1), dydt(nSp + 1), yScale(nSp + 1);

    forAll(rho, c)
    {
        cell_ = c;

        forAll(N_, i)
        {
            y[i] = N_[i].oldTime().internalField()[c];
        }

        y[nSp] =
            1.5*y[eIndex_]*TeF.oldTime().internalField()[c]*kB/e;

        if (debugCharge_)
        {
            // Net charge rate of the chemistry and of the walls
            const bool wallSave = wall_;
            scalarField f(nSp + 1);

            wall_ = false;
            coupledDerivatives(t0, y, f);
            scalar qChem = 0, qPos = 0;
            forAll(N_, i)
            {
                qChem += z(i)*f[i];
                qPos += z(i) > 0 ? z(i)*y[i] : 0;
            }

            wall_ = wallSave;
            coupledDerivatives(t0, y, f);
            scalar qAll = 0, q = 0;
            forAll(N_, i)
            {
                qAll += z(i)*f[i];
                q += z(i)*y[i];
            }

            Info<< "zeroD coupled: t " << t0 << " charge " << q
                << " positive " << qPos << " dq/dt chemistry " << qChem
                << " walls " << qAll - qChem << endl;
        }

        scalar x = t0;
        scalar h = hEst_ > 0 ? min(hEst_, dt) : dt*1e-3;
        label nSteps = 0;

        while (x < t1*(1 - 1e-14) && x < t1 - 1e-300)
        {
            coupledDerivatives(x, y, dydt);

            forAll(y, i)
            {
                yScale[i] = mag(y[i]) + absTol_;
            }

            if (x + h > t1)
            {
                h = t1 - x;
            }

            scalar hDid, hNext;
            rosenbrockStep(x, y, dydt, yScale, h, hDid, hNext);

            h = hNext;

            if (++nSteps > 1000000)
            {
                FatalErrorIn("zeroD::coupledStep")
                    << "more than 1e6 steps in one time step"
                    << exit(FatalError);
            }
        }

        hEst_ = h;

        // Negative densities within the tolerance left by the method
        forAll(y, i)
        {
            if (y[i] < 0)
            {
                if (y[i] < -10*absTol_)
                {
                    WarningIn("zeroD::coupledStep")
                        << "density " << y[i] << " of variable " << i
                        << " at t = " << t1 << endl;
                }

                y[i] = 0;
            }
        }

        // Densities, mass fractions and Te of the end of the step
        TeF.internalField()[c] = coupledTe(t1, y);

        scalar yt = 0;

        forAll(N_, i)
        {
            if (i != bIndex_)
            {
                volScalarField& yi = thermo_.composition().Y(i);

                yi.internalField()[c] = y[i]*W(i)/(Av*rho[c]);
                yt += yi.internalField()[c];
            }
        }

        thermo_.composition().Y(bIndex_).internalField()[c] = 1 - yt;

        // Energy carried to the walls by the electrons at the end of the
        // step, for output
        if (wall_)
        {
            scalarList kIon, kNeutral;
            scalarField n(nSp);

            forAll(n, i)
            {
                n[i] = y[i];
            }

            scalar Ng = NG_.internalField()[c];

            if (Ng <= 0)
            {
                Ng = sum(n);
            }

            wallFrequencies
            (
                n, TeF.internalField()[c], thermo_.T().internalField()[c],
                Ng, kIon, kNeutral
            );

            scalar eRate = 0;

            forAll(wallIon_, k)
            {
                eRate += kIon[k]*n[wallIon_[k]];
            }

            electronWallEnergyLoss_().internalField()[c] =
                2*kB*TeF.internalField()[c]*eRate;
        }
    }

    forAll(N_, i)
    {
        volScalarField& yi = thermo_.composition().Y(i);
        yi.correctBoundaryConditions();
        N_[i] == thermo_.rho()*yi*plasmaConstants::A/W(i);
    }

    TeF.correctBoundaryConditions();
    updateTemperature();
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

    if (coupled_)
    {
        // Species and electron energy over the step as one stiff system
        thermo_.Te().oldTime();

        coupledStep(chemistry, dt);
    }
    else
    {
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
    }

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
