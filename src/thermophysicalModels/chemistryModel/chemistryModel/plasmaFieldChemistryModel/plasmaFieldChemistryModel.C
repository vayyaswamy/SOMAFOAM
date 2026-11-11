/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "plasmaFieldChemistryModel.H"
#include "electronImpactInelasticReaction.H"
#include "electronImpactElasticReaction.H"
#include "TabularReactionRate.H"

// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace
{
    // Factor c^e of a reactant in omega, as ODEChemistryModel::omegaForward:
    // negative concentrations count as zero, and with e < 1 a concentration
    // at or below ROOTVSMALL gives zero
    inline Foam::scalar reactantFactor(const Foam::scalar c, const Foam::scalar e)
    {
        const Foam::scalar cp = Foam::max(c, 0.0);

        if (e == 1.0)
        {
            return cp;
        }
        if (e < 1.0 && cp <= Foam::ROOTVSMALL)
        {
            return 0.0;
        }

        return Foam::pow(cp, e);
    }

    // d(c^e)/dc as in the linearization of ODEChemistryModel::calculateWcf
    inline Foam::scalar reactantFactorDerivative
    (
        const Foam::scalar c,
        const Foam::scalar e
    )
    {
        const Foam::scalar cp = Foam::max(c, 0.0);

        if (e == 1.0)
        {
            return 1.0;
        }
        if (e < 1.0)
        {
            return
                cp > Foam::ROOTVSMALL
              ? e*Foam::pow(cp + Foam::VSMALL, e - 1.0)
              : 0.0;
        }

        return e*Foam::pow(cp, e - 1.0);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::setUp()
{
    const PtrList<Reaction<ThermoType> >& reactions = this->reactions_;
    const label nR = reactions.size();
    const label nS = this->nSpecie_;

    lhsSpecie_.setSize(nR);
    lhsCoeff_.setSize(nR);
    lhsExp_.setSize(nR);
    rhsSpecie_.setSize(nR);
    rhsCoeff_.setSize(nR);
    lhsNet_.setSize(nR);
    tableT_.setSize(nR);
    tableLnK_.setSize(nR);
    fieldRate_.setSize(nR);
    cacheT_.setSize(nR);
    cacheP_.setSize(nR);
    cacheKf_.setSize(nR);
    cacheDkf_.setSize(nR);
    cacheDerivative_.setSize(nR, false);

    reversible_ = false;

    const scalarField cA(nS, 1.0);
    const scalarField cB(nS, 2.0);

    typedef electronImpactInelasticReaction<ThermoType, TabularReactionRate>
        inelasticTabular;
    typedef electronImpactElasticReaction<ThermoType, TabularReactionRate>
        elasticTabular;

    forAll(reactions, i)
    {
        const Reaction<ThermoType>& R = reactions[i];

        const label nL = R.lhs().size();
        const label nRhs = R.rhs().size();

        lhsSpecie_[i].setSize(nL);
        lhsCoeff_[i].setSize(nL);
        lhsExp_[i].setSize(nL);
        lhsNet_[i].setSize(nL);
        rhsSpecie_[i].setSize(nRhs);
        rhsCoeff_[i].setSize(nRhs);

        forAll(R.lhs(), j)
        {
            lhsSpecie_[i][j] = R.lhs()[j].index;
            lhsCoeff_[i][j] = R.lhs()[j].stoichCoeff;
            lhsExp_[i][j] = R.lhs()[j].exponent;
        }

        forAll(R.rhs(), j)
        {
            rhsSpecie_[i][j] = R.rhs()[j].index;
            rhsCoeff_[i][j] = R.rhs()[j].stoichCoeff;
        }

        forAll(R.lhs(), j)
        {
            const label sj = lhsSpecie_[i][j];
            scalar net = 0;

            forAll(R.rhs(), k)
            {
                if (rhsSpecie_[i][k] == sj)
                {
                    net += rhsCoeff_[i][k];
                }
            }
            forAll(R.lhs(), k)
            {
                if (lhsSpecie_[i][k] == sj)
                {
                    net -= lhsCoeff_[i][k];
                }
            }

            lhsNet_[i][j] = net;
        }

        // Reverse rates at gas, ion and electron temperatures
        const scalar Ts[3] = {300.0, 1e4, 5e4};

        for (label t = 0; t < 3; t++)
        {
            const scalar kf = R.kf(Ts[t], 1e2, cA);

            if (R.kr(kf, Ts[t], 1e2, cA) != 0)
            {
                reversible_ = true;
            }
        }

        fieldRate_[i] = R.kf(1e4, 1e2, cA) == R.kf(1e4, 1e2, cB);

        const TabularReactionRate* table = NULL;

        if (const inelasticTabular* r = dynamic_cast<const inelasticTabular*>(&R))
        {
            table = &r->rate();
        }
        else if (const elasticTabular* r = dynamic_cast<const elasticTabular*>(&R))
        {
            table = &r->rate();
        }

        if (table && table->increasing() && table->x().size() > 1)
        {
            tableT_[i] = table->x();
            tableLnK_[i] = table->y();
        }
    }

    if (reversible_)
    {
        Info<< "plasmaFieldChemistryModel: the mechanism has reverse rates;"
            << " the source terms are calculated by ODEChemistryModel"
            << endl;
    }

    setUp_ = true;
}


template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::concentrations
(
    const scalarField& rho
)
{
    const label nS = this->nSpecie_;
    const label nC = rho.size();

    if (c_.size() != nS || kf_.size() != nC)
    {
        c_.setSize(nS);
        dcdt_.setSize(nS);

        for (label s = 0; s < nS; s++)
        {
            c_.set(s, new scalarField(nC));
            dcdt_.set(s, new scalarField(nC));
        }

        kf_.setSize(nC);
        dkfdTe_.setSize(nC);
        prodC_.setSize(nC);
    }

    for (label s = 0; s < nS; s++)
    {
        const scalarField& Y = this->Y_[s].internalField();
        const scalar W = this->specieThermo_[s].W();
        scalarField& c = c_[s];

        forAll(c, celli)
        {
            c[celli] = rho[celli]*Y[celli]/W;
        }
    }
}


template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::forwardRate
(
    const label i,
    const scalarField& Tr,
    const scalarField& p,
    const bool withDerivative
)
{
    // Coefficients of the last evaluation, if the temperatures and pressures
    // are the same (the same values as evaluating them again)
    if
    (
        fieldRate_[i]
     && cacheT_[i].size() == Tr.size()
     && (cacheDerivative_[i] || !withDerivative)
     && cacheT_[i] == Tr
     && cacheP_[i] == p
    )
    {
        kf_ = cacheKf_[i];

        if (cacheDerivative_[i])
        {
            dkfdTe_ = cacheDkf_[i];
        }

        return;
    }

    const Reaction<ThermoType>& R = this->reactions_[i];
    const scalarList& x = tableT_[i];
    const scalarList& y = tableLnK_[i];
    const label n = x.size();

    // Tabulated rates give their derivative in every cell
    bool derivative = withDerivative || n > 1;

    // Concentrations of the cell, for rates that depend on them
    scalarField cCell(this->nSpecie_, 1.0);

    label lo = 0;

    forAll(Tr, celli)
    {
        const scalar T = Tr[celli];

        if (n > 1 && T >= x[0] && T < x[n-1])
        {
            // The segment x[lo] <= T < x[lo + 1] of TabularReactionRate,
            // first tried at the segment of the previous cell
            if (!(T >= x[lo] && T < x[lo + 1]))
            {
                lo = 0;
                label hi = n - 1;

                while (hi - lo > 1)
                {
                    const label mid = (lo + hi)/2;

                    if (T >= x[mid])
                    {
                        lo = mid;
                    }
                    else
                    {
                        hi = mid;
                    }
                }
            }

            const label hi = lo + 1;

            if (TabularReactionRate::linearInterpolation)
            {
                // k linearly in Te, as TabularReactionRate
                const scalar klo = exp(y[lo]);
                const scalar dk = (exp(y[hi]) - klo)/(x[hi] - x[lo]);

                kf_[celli] = klo + (exp(y[hi]) - klo)*(T - x[lo])/(x[hi] - x[lo]);
                dkfdTe_[celli] = dk;

                continue;
            }

            const scalar slope = (y[hi] - y[lo])/(x[hi] - x[lo]);

            // The interpolation of TabularReactionRate, term by term
            const scalar k =
                exp
                (
                    y[lo]
                  + (y[hi] - y[lo])
                   *(T - x[lo])
                   /(x[hi] - x[lo])
                );

            kf_[celli] = k;
            dkfdTe_[celli] = k*slope;

            continue;
        }

        if (!fieldRate_[i])
        {
            for (label s = 0; s < this->nSpecie_; s++)
            {
                cCell[s] = c_[s][celli];
            }
        }

        kf_[celli] = R.kf(T, p[celli], cCell);

        if (withDerivative)
        {
            dkfdTe_[celli] =
                (
                    R.kf(T + 1e-08, p[celli], cCell)
                  - R.kf(T - 1e-08, p[celli], cCell)
                )/2e-08;
        }
        else
        {
            // Outside the table without the derivative
            derivative = false;
        }
    }

    if (fieldRate_[i])
    {
        cacheT_[i] = Tr;
        cacheP_[i] = p;
        cacheKf_[i] = kf_;
        cacheDerivative_[i] = derivative;

        if (derivative)
        {
            cacheDkf_[i] = dkfdTe_;
        }
    }
}


template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::reactantProduct
(
    const label i
)
{
    const labelList& sp = lhsSpecie_[i];
    const scalarList& ex = lhsExp_[i];

    prodC_ = 1.0;

    forAll(sp, j)
    {
        const scalarField& c = c_[sp[j]];
        const scalar e = ex[j];

        forAll(prodC_, celli)
        {
            prodC_[celli] *= reactantFactor(c[celli], e);
        }
    }
}


template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::addElectronTerms
(
    const label i,
    const label kind,
    const scalarField& rho,
    const scalarField& T,
    const scalarField& Te
)
{
    // Molar mass of the electrons [kg/kmol], as ODEChemistryModel
    const scalar me = this->electronW();

    if (kind == ODEChemistryModel<CompType, ThermoType>::electronInelastic)
    {
        const scalar slE = lhsCoeff_[i][0]*this->reactions_[i].deltaE()*me;

        forAll(prodC_, celli)
        {
            this->eChemSource_[celli] += slE*kf_[celli]*prodC_[celli];
            this->dEChemSourceDTe_[celli] += slE*dkfdTe_[celli]*prodC_[celli];
        }
    }
    else
    {
        const label sie = lhsSpecie_[i][0];
        const label si = lhsSpecie_[i][1];
        const scalar f =
            3.0*me*me*1.38064852E-23/this->specieThermo_[si].W();
        const scalarField& Ye = this->Y_[sie].internalField();
        scalarField& nu = this->collFreq_[sie];

        forAll(prodC_, celli)
        {
            const scalar omega = kf_[celli]*prodC_[celli];
            const scalar dOmega = dkfdTe_[celli]*prodC_[celli];
            const scalar dT = Te[celli] - T[celli];

            this->eChemSource_[celli] += f*omega*dT;
            nu[celli] += omega*me/Ye[celli]/rho[celli];
            this->dEChemSourceDTe_[celli] += f*(dOmega*dT + omega);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::plasmaFieldChemistryModel<CompType, ThermoType>::plasmaFieldChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    ODEChemistryModel<CompType, ThermoType>
    (
        mesh,
        obj,
        compTypeName,
        thermoTypeName
    ),
    setUp_(false),
    reversible_(false)
{
    Info<< "plasmaFieldChemistryModel: number of species = "
        << this->nSpecie_ << " and reactions = " << this->nReaction_
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::plasmaFieldChemistryModel<CompType, ThermoType>::
~plasmaFieldChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::calculateWcf()
{
    if (this->chemistry_ && !setUp_)
    {
        setUp();
    }

    if (!this->chemistry_ || reversible_)
    {
        ODEChemistryModel<CompType, ThermoType>::calculateWcf();
        return;
    }

    typedef ODEChemistryModel<CompType, ThermoType> base;

    const scalarField rho(this->thermo().rho()().internalField());
    const label nC = rho.size();
    const label nS = this->nSpecie_;

    for (label s = 0; s < nS; s++)
    {
        this->RR_[s].setSize(nC);
        this->dRRDi_[s].setSize(nC);
        this->collFreq_[s].setSize(nC);

        this->dRRDi_[s] = 0.0;
        this->collFreq_[s] = 0.0;
    }

    this->eChemSource_.setSize(nC);
    this->dEChemSourceDTe_.setSize(nC);
    this->eChemSource_ = 0.0;
    this->dEChemSourceDTe_ = 0.0;

    concentrations(rho);

    for (label s = 0; s < nS; s++)
    {
        dcdt_[s] = 0.0;
    }

    const scalarField& T = this->thermo().T().internalField();
    const scalarField& p = this->thermo().p().internalField();
    const scalarField& Te = this->thermo().Te().internalField();
    const scalarField& Tion = this->thermo().Tion().internalField();

    const labelList& kinds = this->reactionKinds();

    forAll(this->reactions_, i)
    {
        const label kind = kinds[i];

        const bool electron =
            kind == base::electronInelastic || kind == base::electronElastic;

        forwardRate
        (
            i,
            electron ? Te : kind == base::ionIrreversible ? Tion : T,
            p,
            electron
        );

        reactantProduct(i);

        // Production of the species
        forAll(lhsSpecie_[i], j)
        {
            scalarField& d = dcdt_[lhsSpecie_[i][j]];
            const scalar sl = lhsCoeff_[i][j];

            forAll(d, celli)
            {
                d[celli] -= sl*kf_[celli]*prodC_[celli];
            }
        }

        forAll(rhsSpecie_[i], j)
        {
            scalarField& d = dcdt_[rhsSpecie_[i][j]];
            const scalar sr = rhsCoeff_[i][j];

            forAll(d, celli)
            {
                d[celli] += sr*kf_[celli]*prodC_[celli];
            }
        }

        // Linearization: derivative of the production of each reactant in
        // its own concentration
        const label nL = lhsSpecie_[i].size();

        forAll(lhsSpecie_[i], j)
        {
            const scalar net = lhsNet_[i][j];

            if (net == 0)
            {
                continue;
            }

            scalarField& dR = this->dRRDi_[lhsSpecie_[i][j]];

            forAll(dR, celli)
            {
                scalar d = kf_[celli];

                for (label k = 0; k < nL; k++)
                {
                    const scalar ck = c_[lhsSpecie_[i][k]][celli];

                    d *=
                        k == j
                      ? reactantFactorDerivative(ck, lhsExp_[i][k])
                      : reactantFactor(ck, lhsExp_[i][k]);
                }

                dR[celli] += net*d;
            }
        }

        if (electron)
        {
            addElectronTerms(i, kind, rho, T, Te);
        }
        else if (kind == base::ionElastic)
        {
            const label sii = lhsSpecie_[i][0];
            const scalar W = this->specieThermo_[sii].W();
            const scalarField& Yi = this->Y_[sii].internalField();
            scalarField& nu = this->collFreq_[sii];

            forAll(nu, celli)
            {
                nu[celli] += kf_[celli]*prodC_[celli]*W/Yi[celli]/rho[celli];
            }
        }
    }

    for (label s = 0; s < nS; s++)
    {
        this->RR_[s] = dcdt_[s]*this->specieThermo_[s].W();
    }
}


template<class CompType, class ThermoType>
void Foam::plasmaFieldChemistryModel<CompType, ThermoType>::
calculateElectronEnergySource()
{
    const label nC = this->mesh().nCells();

    if
    (
        !this->chemistry_
     || !setUp_
     || reversible_
     || this->eChemSource_.size() != nC
     || kf_.size() != nC
    )
    {
        calculateWcf();
        return;
    }

    typedef ODEChemistryModel<CompType, ThermoType> base;

    const scalarField rho(this->thermo().rho()().internalField());

    const labelList& kinds = this->reactionKinds();

    this->eChemSource_ = 0.0;
    this->dEChemSourceDTe_ = 0.0;

    forAll(this->reactions_, i)
    {
        if (kinds[i] == base::electronElastic)
        {
            this->collFreq_[lhsSpecie_[i][0]] = 0.0;
        }
    }

    concentrations(rho);

    const scalarField& T = this->thermo().T().internalField();
    const scalarField& p = this->thermo().p().internalField();
    const scalarField& Te = this->thermo().Te().internalField();

    forAll(this->reactions_, i)
    {
        const label kind = kinds[i];

        if (kind != base::electronInelastic && kind != base::electronElastic)
        {
            continue;
        }

        forwardRate(i, Te, p, true);
        reactantProduct(i);
        addElectronTerms(i, kind, rho, T, Te);
    }
}


// ************************************************************************* //
