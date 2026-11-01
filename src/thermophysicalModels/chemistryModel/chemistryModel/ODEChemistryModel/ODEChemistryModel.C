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

#include "ODEChemistryModel.H"
#include "chemistrySolver.H"
#include "reactingMixture.H"
#include "TabularReactionRate.H"


// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace
{
    // pow(x, e) without the library call for the stoichiometric exponents
    // 1 and 0, where it is exact: pow(x, 1) = x and pow(x, 0) = 1 for any x
    inline Foam::scalar powE(const Foam::scalar x, const Foam::scalar e)
    {
        if (e == 1.0)
        {
            return x;
        }
        if (e == 0.0)
        {
            return 1.0;
        }
        return Foam::pow(x, e);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ODEChemistryModel<CompType, ThermoType>::ODEChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    CompType(mesh, obj, thermoTypeName),

    ODE(),

    Y_(this->thermo().composition().Y()),

    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),

    solver_
    (
        chemistrySolver<CompType, ThermoType>::New
        (
            *this,
            compTypeName,
            thermoTypeName
        )
    ),

    RR_(nSpecie_),
    coeffs_(nSpecie_ + 4),
	eChemSource_(mesh.nCells(), 0.0),
    dEChemSourceDTe_(mesh.nCells(), 0.0),
	collFreq_(nSpecie_),
    dRRDi_(nSpecie_),
    rateCache_(),
    rateCacheCells_(-1),
    cacheableRate_()
{
    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );
        dRRDi_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );
    }
    forAll(collFreq_, fieldI)
    {
        collFreq_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );
	}
    // Interpolation of the tabulated rates (TabularReactionRate) between
    // the rows of a table: rateInterpolation logarithmic (default, ln k
    // linearly in Te) or linear (k linearly in Te, and so in the mean
    // energy 3/2 kB Te)
    const word rateInterpolation
    (
        this->template lookupOrDefault<word>("rateInterpolation", "logarithmic")
    );

    if (rateInterpolation != "logarithmic" && rateInterpolation != "linear")
    {
        FatalErrorIn("ODEChemistryModel::ODEChemistryModel")
            << "rateInterpolation " << rateInterpolation
            << " in chemistryProperties: logarithmic or linear"
            << exit(FatalError);
    }

    TabularReactionRate::linearInterpolation = rateInterpolation == "linear";

    Info<< "ODEChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_
        << ", tabulated rates interpolated " << rateInterpolation << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ODEChemistryModel<CompType, ThermoType>::~ODEChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalarField Foam::ODEChemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    const scalar Te,
    const scalar Tion
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalarField om(nEqns(), 0.0);

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, Te, Tion, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return om;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    const scalar Te,
    const scalar Tion,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    return omegaKind
    (
        kindOf(R.type()), R, c, T, p, Te, Tion, pf, cf, lRef, pr, cr, rRef
    );
}


template<class CompType, class ThermoType>
Foam::label Foam::ODEChemistryModel<CompType, ThermoType>::kindOf
(
    const word& type
)
{
    if
    (
        type == "electronImpactInelasticArrheniusReaction"
     || type == "electronImpactInelasticTabularReaction"
    )
    {
        return electronInelastic;
    }
    if (type == "ionIrreversibleArrheniusReaction")
    {
        return ionIrreversible;
    }
    if
    (
        type == "electronImpactElasticArrheniusReaction"
     || type == "electronImpactElasticTabularReaction"
    )
    {
        return electronElastic;
    }
    if (type == "ionElasticArrheniusReaction")
    {
        return ionElastic;
    }
    return otherReaction;
}


template<class CompType, class ThermoType>
const Foam::labelList&
Foam::ODEChemistryModel<CompType, ThermoType>::reactionKinds() const
{
    if (reactionKind_.size() != reactions_.size())
    {
        reactionKind_.setSize(reactions_.size());
        forAll(reactions_, i)
        {
            reactionKind_[i] = kindOf(reactions_[i].type());
        }
    }
    return reactionKind_;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::omegaKind
(
    const label kind,
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    const scalar Te,
    const scalar Tion,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{

    const scalar Tr = rateTemperature(kind, T, Te, Tion);

    const scalar kf = R.kf(Tr, p, c);
    const scalar kr = R.kr(kf, Tr, p, c);

    return omegaFromRates(R, c, kf, kr, pf, cf, lRef, pr, cr, rRef);
}


template<class CompType, class ThermoType>
Foam::scalar
Foam::ODEChemistryModel<CompType, ThermoType>::rateTemperature
(
    const label kind,
    const scalar T,
    const scalar Te,
    const scalar Tion
)
{
    // Electron-impact rates, elastic and inelastic, are functions of the
    // electron temperature
    if (kind == electronInelastic || kind == electronElastic)
    {
        return Te;
    }
    else if (kind == ionIrreversible)
    {
        return Tion;
    }

    return T;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::cachedRate
(
    const label i,
    const label celli,
    const label slot,
    const scalar T,
    const scalar p,
    const scalarField& c,
    scalar& kf,
    scalar& kr
) const
{
    const Reaction<ThermoType>& R = reactions_[i];

    if (!cacheableRate_[i])
    {
        kf = R.kf(T, p, c);
        kr = R.kr(kf, T, p, c);
        return;
    }

    scalar* e = &rateCache_[4*((celli*nReaction_ + i)*rateSlots + slot)];

    if (e[0] == T && e[1] == p)
    {
        kf = e[2];
        kr = e[3];
        return;
    }

    kf = R.kf(T, p, c);
    kr = R.kr(kf, T, p, c);

    e[0] = T;
    e[1] = p;
    e[2] = kf;
    e[3] = kr;
}


template<class CompType, class ThermoType>
Foam::scalar
Foam::ODEChemistryModel<CompType, ThermoType>::omegaForward
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar kf
) const
{
    // The forward part of omegaFromRates, term by term
    label Nl = R.lhs().size();

    label slRef = 0;
    label lRef = R.lhs()[slRef].index;

    scalar pf = kf;
    for (label s=1; s<Nl; s++)
    {
        label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            scalar exp = R.lhs()[slRef].exponent;
            pf *= powE(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= powE(max(0.0, c[si]), exp);
        }
    }
    const scalar cf = max(0.0, c[lRef]);

    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > ROOTVSMALL)
            {
                pf *= powE(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= powE(cf, exp - 1.0);
        }
    }

    return pf*cf;
}


template<class CompType, class ThermoType>
Foam::scalar
Foam::ODEChemistryModel<CompType, ThermoType>::omegaFromRates
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar kf,
    const scalar kr,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    pf = 1.0;
    pr = 1.0;

    label Nl = R.lhs().size();
    label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            scalar exp = R.lhs()[slRef].exponent;
            pf *= powE(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= powE(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);


    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > ROOTVSMALL)
            {
                pf *= powE(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= powE(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            scalar exp = R.rhs()[srRef].exponent;
            pr *= powE(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            scalar exp = R.rhs()[s].exponent;
            pr *= powE(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        scalar exp = R.rhs()[srRef].exponent;
        if (exp<1.0)
        {
            if (cr>ROOTVSMALL)
            {
                pr *= powE(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= powE(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    scalar T = c[nSpecie_];
    scalar p = c[nSpecie_ + 1];
    scalar Te = c[nSpecie_ + 2];
    scalar Tion = c[nSpecie_ + 3];

    dcdt = omega(c, T, p, Te, Tion);

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        scalar W = specieThermo_[i].W();
        cSum += c[i];
        rho += W*c[i];
    }
    scalar mw = rho/cSum;
    scalar cp = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        scalar cpi = specieThermo_[i].cp(T);
        scalar Xi = c[i]/rho;
        cp += Xi*cpi;
    }
    cp /= mw;

    scalar dT = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        scalar hi = specieThermo_[i].h(T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    // limit the time-derivative, this is more stable for the ODE
    // solver when calculating the allowed time step
    scalar dtMag = min(500.0, mag(dT));
    dcdt[nSpecie_] = -dT*dtMag/(mag(dT) + SMALL);

    // dp/dt = ...
    dcdt[nSpecie_+1] = 0.0;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    scalar T = c[nSpecie_];
    scalar p = c[nSpecie_ + 1];
    scalar Te = c[nSpecie_ + 2];
    scalar Tion = c[nSpecie_ + 3];

    scalarField c2(nSpecie_, 0.0);
    for (label i=0; i<nSpecie_; i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    dcdt = omega(c2, T, p, Te, Tion);

    for (label ri=0; ri<reactions_.size(); ri++)
    {
        const Reaction<ThermoType>& R = reactions_[ri];

    	scalar kf0, kr0;
		if ((R.type() == "electronImpactInelasticArrheniusReaction") || (R.type() == "electronImpactInelasticTabularReaction"))
		{
		  kf0 = R.kf(Te, p, c2);
		  kr0 = R.kr(kf0, Te, p, c2);
		}
		else if ((R.type() == "ionIrreversibleArrheniusReaction"))
		{
		  kf0 = R.kf(Tion, p, c2);
		  kr0 = R.kr(kf0, Tion, p, c2);
		}
		else
		{
		  kf0 = R.kf(T, p, c2);
		  kr0 = R.kr(kf0, T, p, c2);
		}

        forAll(R.lhs(), j)
        {
            label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c2[si]>ROOTVSMALL)
                        {
                            kf *= el*powE(c2[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*powE(c2[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= powE(c2[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar er = R.rhs()[i].exponent;
                if (i==j)
                {
                    if (er<1.0)
                    {
                        if (c2[si]>ROOTVSMALL)
                        {
                            kr *= er*powE(c2[si] + VSMALL, er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*powE(c2[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= powE(c2[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] -= sr*kr;
            }
        }


		// calculate the dcdT elements numerically
		scalar delta = 1.0e-8;
		scalarField dcdT0, dcdT1;
		if ((R.type() == "electronImpactInelasticArrheniusReaction") || (R.type() == "electronImpactInelasticTabularReaction"))
		{
		  dcdT0 = omega(c2, T, p, Te - delta, Tion);
		  dcdT1 = omega(c2, T, p, Te + delta, Tion);
		}
		else if ((R.type() == "ionIrreversibleArrheniusReaction"))
		{
		  dcdT0 = omega(c2, T, p, Te, Tion - delta);
		  dcdT1 = omega(c2, T, p, Te, Tion + delta);
		}
		else
		{
		  dcdT0 = omega(c2, T - delta, p, Te, Tion);
		  dcdT1 = omega(c2, T + delta, p, Te, Tion);
		}
	//    scalarField dcdT0 = omega(c2, T - delta, p, Te, Tion);
	//    scalarField dcdT1 = omega(c2, T + delta, p, Te, Tion);

		for (label ni=0; ni<nEqns(); ni++)
		{
		    dfdc[ni][nSpecie()] = 0.5*(dcdT1[ni] - dcdT0[ni])/delta;
		}
    }
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::tc() const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc();

    label nReaction = reactions_.size();

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];
            scalar Tei = this->thermo().Te()[celli];
            scalar Tioni = this->thermo().Tion()[celli];
            scalarField c(nSpecie_);
            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                cSum += c[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];

                omega
                (
                    R, c, Ti, pi, Tei, Tioni, pf, cf, lRef, pr, cr, rRef
                );

                forAll(R.rhs(), s)
                {
                    scalar sr = R.rhs()[s].stoichCoeff;
                    tc[celli] += sr*pf*cf;
                }
            }
            tc[celli] = nReaction*cSum/tc[celli];
        }
    }


    ttc().correctBoundaryConditions();

    return ttc;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                scalar hi = specieThermo_[i].Hc();
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }

    return tSh;
}

template<class CompType, class ThermoType>
bool Foam::ODEChemistryModel<CompType, ThermoType>::electronElasticLoss() const
{
    forAll(reactions_, i)
    {
        const word& type = reactions_[i].type();

        if
        (
            type == "electronImpactElasticArrheniusReaction"
         || type == "electronImpactElasticTabularReaction"
        )
        {
            return true;
        }
    }

    return false;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::eChemSourceExplicit() const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    tmp<volScalarField> tShe
    (
        new volScalarField
        (
            IOobject
            (
                "She",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy*dimMass/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& She = tShe();

        const scalarField& T = this->thermo().T();
        const scalarField& p = this->thermo().p();
        const scalarField& Te = this->thermo().Te();
        const scalarField& Tion = this->thermo().Tion();

        forAll(She, celli)
        {
	    	const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = p[celli];
            const scalar Tei = Te[celli];
            const scalar Tioni = Tion[celli];

            scalarField c(nSpecie_, 0.0);
            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
            }
        
            forAll(reactions_, i)
            {
	        	const Reaction<ThermoType>& R = reactions_[i];

                scalar omegai = omega
                (
	       			R, c, Ti, pi, Tei, Tioni, pf, cf, lRef, pr, cr, rRef
                );
				if ((R.type() == "electronImpactInelasticArrheniusReaction") || (R.type() == "electronImpactInelasticTabularReaction"))
				{
		        	She[celli] += omegai*R.deltaE()*5.48579E-4; //*dimensionedScalar("one",dimEnergy,1.0); // in J.kg/m^3/s, electron weight hard coded 
				}
				if ((R.type() == "electronImpactElasticArrheniusReaction") || (R.type() == "electronImpactElasticTabularReaction"))
				{
				  // getting index of neutral species
					const label si = R.lhs()[1].index; 

					She[celli] += 3.0*omegai*(Tei-Ti)*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W(); //electron weight square
				  // need to multiply by electron mass and divide by other species mass
				}
            }              
        }
    }
    tShe().correctBoundaryConditions();  

    return tShe;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class ThermoType>
Foam::label Foam::ODEChemistryModel<CompType, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 4;
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::calculate()
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    for (label i=0; i<nSpecie_; i++)
    {
        RR_[i].setSize(rho.size());
        dRRDi_[i].setSize(rho.size());
        collFreq_[i].setSize(rho.size());
    }

    // The number of cells changes with adaptive mesh refinement
    eChemSource_.setSize(rho.size());
    dEChemSourceDTe_.setSize(rho.size());

    Info << "In calculate" << endl;

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0.0;
                dRRDi_[i][celli] = 0.0;
            }

			eChemSource_[celli] = 0.0;

            dEChemSourceDTe_[celli] = 0.0;

            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];
            scalar Tei = this->thermo().Te()[celli];
            scalar Tioni = this->thermo().Tion()[celli];

            scalarField c(nSpecie_);
            scalarField c2(nSpecie_, 0.0);
            scalarField dcdt(nEqns(), 0.0);

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                c2[i] = max(c[i], 0.0);
            }

			scalar pf, cf, pr, cr;
			label lRef, rRef;

			forAll(reactions_, i)
			{
				const Reaction<ThermoType>& R = reactions_[i];

				scalar omegai = omega
				(
				    R, c, Ti, pi, Tei, Tioni, pf, cf, lRef, pr, cr, rRef
				);

                scalar omegaiplus = omega
                (
                    R, c, Ti, pi, Tei + 1e-08, Tioni, pf, cf, lRef, pr, cr, rRef
                );

                scalar omegaiminus = omega
                (
                    R, c, Ti, pi, Tei - 1e-08, Tioni, pf, cf, lRef, pr, cr, rRef
                );

                scalar dOmegaiDTe = (omegaiplus - omegaiminus) / 2e-08; 

                scalar kf0, kr0;


                if ((R.type() == "electronImpactInelasticArrheniusReaction") || (R.type() == "electronImpactInelasticTabularReaction"))
                {
                    kf0 = R.kf(Tei, pi, c2);
                    kr0 = R.kr(kf0, Tei, pi, c2);
                }
                else if ((R.type() == "ionIrreversibleArrheniusReaction"))
                {
                    kf0 = R.kf(Tioni, pi, c2);
                    kr0 = R.kr(kf0, Tioni, pi, c2);
                }
                else
                {
                    kf0 = R.kf(Ti, pi, c2);
                    kr0 = R.kr(kf0, Ti, pi, c2);
                }

                forAll(R.lhs(), j)
                {
                    label sj = R.lhs()[j].index;
                    scalar kf = kf0;
                    forAll(R.lhs(), i)
                    {
                        label si = R.lhs()[i].index;
                        scalar el = R.lhs()[i].exponent;
                        if (i == j)
                        {
                            if (el < 1.0)
                            {
                                if (c2[si]>ROOTVSMALL)
                                {
                                    kf *= el*powE(c2[si] + VSMALL, el - 1.0);
                                }
                                else
                                {
                                    kf = 0.0;
                                }
                            }
                            else
                            {
                                kf *= el*powE(c2[si], el - 1.0);
                            }
                        }
                        else
                        {
                            kf *= powE(c2[si], el);
                        }
                    }

                    forAll(R.lhs(), i)
                    {
                        label si = R.lhs()[i].index;
                        scalar sl = R.lhs()[i].stoichCoeff;
                        
                        if (si == sj)
                        {
                            dRRDi_[si][celli] -= sl*kf;
                        }
                    }
                    forAll(R.rhs(), i)
                    {
                        label si = R.rhs()[i].index;
                        scalar sr = R.rhs()[i].stoichCoeff;

                        if (si == sj)
                        {
                            dRRDi_[si][celli] += sr*kf;
                        }
                    }
                }

                forAll(R.rhs(), j)
                {
                    label sj = R.rhs()[j].index;
                    scalar kr = kr0;
                    forAll(R.rhs(), i)
                    {
                        label si = R.rhs()[i].index;
                        scalar er = R.rhs()[i].exponent;
                        if (i==j)
                        {
                            if (er<1.0)
                            {
                                if (c2[si]>ROOTVSMALL)
                                {
                                    kr *= er*powE(c2[si] + VSMALL, er - 1.0);
                                }
                                else
                                {
                                    kr = 0.0;
                                }
                            }
                            else
                            {
                                kr *= er*powE(c2[si], er - 1.0);
                            }
                        }
                        else
                        {
                            kr *= powE(c2[si], er);
                        }   
                    }

                    forAll(R.lhs(), i)
                    {
                        label si = R.lhs()[i].index;
                        scalar sl = R.lhs()[i].stoichCoeff;
                        if (si == sj)
                        {
                            dRRDi_[si][celli] += sl*kr;    
                        }
                        
                    }
                    forAll(R.rhs(), i)
                    {
                        label si = R.rhs()[i].index;
                        scalar sr = R.rhs()[i].stoichCoeff;
                        if (si == sj)
                        {
                            dRRDi_[si][celli] -= sr*kr;    
                        }
                    }
                }





				if ((R.type() == "electronImpactInelasticArrheniusReaction") || (R.type() == "electronImpactInelasticTabularReaction"))
				{
				    scalar sl = R.lhs()[0].stoichCoeff;
			    	eChemSource_[celli] += sl*omegai*R.deltaE()*5.48579E-4; //*dimensionedScalar("one",dimEnergy,1.0); // in J.kg/m^3/s, electron weight hard coded 
				    dEChemSourceDTe_[celli] += sl*dOmegaiDTe*R.deltaE()*5.48579E-4;
                }
				if ((R.type() == "electronImpactElasticArrheniusReaction") || (R.type() == "electronImpactElasticTabularReaction"))
				{
				  // getting index of neutral species
					const label si = R.lhs()[1].index; 

					eChemSource_[celli] += 3.0*omegai*(Tei-Ti)*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W(); //electron weight square
				  // need to multiply by electron mass and divide by other species mass
                    dEChemSourceDTe_[celli] += 3.0*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W()*(dOmegaiDTe*(Tei-Ti) + omegai);
				}

				forAll(R.lhs(), s)
				{
				    label si = R.lhs()[s].index;
				    scalar sl = R.lhs()[s].stoichCoeff;
				    dcdt[si] -= sl*omegai;
				}

				forAll(R.rhs(), s)
				{
				    label si = R.rhs()[s].index;
				    scalar sr = R.rhs()[s].stoichCoeff;
				    dcdt[si] += sr*omegai;
				}
			}

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
                Info << "dRR = " << dRRDi_[i][celli] << endl;
            }
        }
    }
}

template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::calculateWcf()
{
    //Info << "Inside calculateWcf" << endl;
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    for (label i=0; i<nSpecie_; i++)
    {
        RR_[i].setSize(rho.size());
        dRRDi_[i].setSize(rho.size());
        collFreq_[i].setSize(rho.size());
    }

    // The number of cells changes with adaptive mesh refinement
    eChemSource_.setSize(rho.size());
    dEChemSourceDTe_.setSize(rho.size());

    if (this->chemistry_)
    {
        const labelList& kinds = reactionKinds();

        // Rates that do not depend on the concentrations can be cached:
        // kf and kr at two concentration vectors agree
        if (cacheableRate_.size() != nReaction_)
        {
            cacheableRate_.setSize(nReaction_);

            const scalarField cA(nSpecie_, 1.0);
            const scalarField cB(nSpecie_, 2.0);

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];
                const scalar T0 = 1e4;
                const scalar p0 = 1e2;
                const scalar kfA = R.kf(T0, p0, cA);
                const scalar kfB = R.kf(T0, p0, cB);

                cacheableRate_[i] =
                    kfA == kfB
                 && R.kr(kfA, T0, p0, cA) == R.kr(kfB, T0, p0, cB);
            }
        }

        if (rateCacheCells_ != rho.size())
        {
            rateCache_.setSize(4*rho.size()*nReaction_*rateSlots);
            rateCache_ = -GREAT;
            rateCacheCells_ = rho.size();
        }

        const scalarField& TCells = this->thermo().T();
        const scalarField& pCells = this->thermo().p();
        const scalarField& TeCells = this->thermo().Te();
        const scalarField& TionCells = this->thermo().Tion();

        // Work arrays, set in every cell
        scalarField c(nSpecie_);
        scalarField dcdt(nEqns(), 0.0);
        scalarField c2(nSpecie_, 0.0);

        forAll(rho, celli)
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0.0;
				collFreq_[i][celli] = 0.0;
                dRRDi_[i][celli] = 0.0;
            }

			eChemSource_[celli] = 0.0;

            dEChemSourceDTe_[celli] = 0.0;

            scalar rhoi = rho[celli];
            scalar Ti = TCells[celli];
            scalar pi = pCells[celli];
            scalar Tei = TeCells[celli];
            scalar Tioni = TionCells[celli];

            dcdt = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                c2[i] = max(c[i], 0.0);
            }

			scalar pf, cf, pr, cr;
			label lRef, rRef;

            

			forAll(reactions_, i)
			{
				const Reaction<ThermoType>& R = reactions_[i];
				const label kind = kinds[i];

				// Rate coefficients from the cache (unchanged when the
				// temperature of the reaction is), as omegaKind
				const scalar Tr = rateTemperature(kind, Ti, Tei, Tioni);

				scalar kfr, krr;
				cachedRate(i, celli, 0, Tr, pi, c, kfr, krr);

				// With kr = 0 the reverse part pr*cr of omega is +0, and
				// pf*cf - 0 = pf*cf
				scalar omegai =
				    krr == 0.0
				  ? omegaForward(R, c, kfr)
				  : omegaFromRates(R, c, kfr, krr, pf, cf, lRef, pr, cr, rRef);

                // Central difference in Te; only the electron-impact
                // rates depend on Te (omegaKind), for the other reactions
                // the difference is exactly zero
                scalar dOmegaiDTe = 0.0;

                if (kind == electronInelastic || kind == electronElastic)
                {
                    scalar kfp, krp, kfm, krm;
                    cachedRate(i, celli, 1, Tei + 1e-08, pi, c, kfp, krp);
                    cachedRate(i, celli, 2, Tei - 1e-08, pi, c, kfm, krm);

                    scalar omegaiplus =
                        krp == 0.0
                      ? omegaForward(R, c, kfp)
                      : omegaFromRates(R, c, kfp, krp, pf, cf, lRef, pr, cr, rRef);

                    scalar omegaiminus =
                        krm == 0.0
                      ? omegaForward(R, c, kfm)
                      : omegaFromRates(R, c, kfm, krm, pf, cf, lRef, pr, cr, rRef);

                    dOmegaiDTe = (omegaiplus - omegaiminus) / 2e-08;
                }

                // Linearization: inelastic at Te, ion reactions at Tion,
                // the others (elastic included) at T
                scalar kf0, kr0;

                const scalar T0 =
                    kind == electronInelastic ? Tei
                  : kind == ionIrreversible ? Tioni
                  : Ti;

                if (cacheableRate_[i] && T0 == Tr)
                {
                    kf0 = kfr;
                    kr0 = krr;
                }
                else
                {
                    cachedRate(i, celli, 3, T0, pi, c2, kf0, kr0);
                }

                forAll(R.lhs(), j)
                {
                    label sj = R.lhs()[j].index;
                    scalar kf = kf0;
                    forAll(R.lhs(), i)
                    {
                        label si = R.lhs()[i].index;
                        scalar el = R.lhs()[i].exponent;
                        // Differentiate the j-th factor only: a species may appear
                        // more than once on the left-hand side (A + A), as in
                        // calculate() and jacobian()
                        if (i == j)
                        {
                            if (el < 1.0)
                            {
                                if (c2[si]>ROOTVSMALL)
                                {
                                    kf *= el*powE(c2[si] + VSMALL, el - 1.0);
                                }
                                else
                                {
                                    kf = 0.0;
                                }
                            }
                            else
                            {
                                kf *= el*powE(c2[si], el - 1.0);
                            }
                        }
                        else
                        {
                            kf *= powE(c2[si], el);
                        }
                    }

                    forAll(R.lhs(), i)
                    {
                        label si = R.lhs()[i].index;
                        scalar sl = R.lhs()[i].stoichCoeff;
                        
                        if (si == sj)
                        {
                            dRRDi_[si][celli] -= sl*kf;
                        }
                    }
                    forAll(R.rhs(), i)
                    {
                        label si = R.rhs()[i].index;
                        scalar sr = R.rhs()[i].stoichCoeff;

                        if (si == sj)
                        {
                            dRRDi_[si][celli] += sr*kf;
                        }
                    }
                    //Info << "dRRDi_[sj][celli]" << dRRDi_[sj][celli] << endl;
                }

                // With kr0 = 0 every reverse term below adds or subtracts +0
                // to entries that are never -0: nothing changes
                forAll(R.rhs(), j)
                {
                    if (kr0 == 0.0)
                    {
                        break;
                    }

                    label sj = R.rhs()[j].index;
                    scalar kr = kr0;
                    forAll(R.rhs(), i)
                    {
                        label si = R.rhs()[i].index;
                        scalar er = R.rhs()[i].exponent;
                        if (i==j)
                        {
                            if (er<1.0)
                            {
                                if (c2[si]>ROOTVSMALL)
                                {
                                    kr *= er*powE(c2[si] + VSMALL, er - 1.0);
                                }
                                else
                                {
                                    kr = 0.0;
                                }
                            }
                            else
                            {
                                kr *= er*powE(c2[si], er - 1.0);
                            }
                        }
                        else
                        {
                            kr *= powE(c2[si], er);
                        }   
                    }

                    forAll(R.lhs(), i)
                    {
                        label si = R.lhs()[i].index;
                        scalar sl = R.lhs()[i].stoichCoeff;
                        if (si == sj)
                        {
                            dRRDi_[si][celli] += sl*kr;    
                        }
                        
                    }
                    forAll(R.rhs(), i)
                    {
                        label si = R.rhs()[i].index;
                        scalar sr = R.rhs()[i].stoichCoeff;
                        if (si == sj)
                        {
                            dRRDi_[si][celli] -= sr*kr;    
                        }
                    }
                }

				if (kind == electronInelastic)
				{
				    scalar sl = R.lhs()[0].stoichCoeff;
			    	eChemSource_[celli] += sl*omegai*R.deltaE()*5.48579E-4; //*dimensionedScalar("one",dimEnergy,1.0); // in J.kg/m^3/s, electron weight hard coded 
				    dEChemSourceDTe_[celli] += sl*dOmegaiDTe*R.deltaE()*5.48579E-4;
                }
				if (kind == electronElastic)
				{
				  // getting index of neutral species
					const label si = R.lhs()[1].index; 
					const label sie = R.lhs()[0].index; 
					eChemSource_[celli] += 3.0*omegai*(Tei-Ti)*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W(); //electron weight square
					collFreq_[sie][celli] += omegai*5.48579E-4/Y_[sie][celli]/rhoi;
                    dEChemSourceDTe_[celli] += 3.0*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W()*(dOmegaiDTe*(Tei-Ti) + omegai);
				  // need to multiply by electron mass and divide by other species mass
				}
				if (kind == ionElastic)
				{
				  // getting index of neutral species
					const label sii = R.lhs()[0].index; 
					collFreq_[sii][celli] += omegai*specieThermo_[sii].W()/Y_[sii][celli]/rhoi;
				  // need to multiply by electron mass and divide by other species mass
				}

				forAll(R.lhs(), s)
				{
				    label si = R.lhs()[s].index;
				    scalar sl = R.lhs()[s].stoichCoeff;
				    dcdt[si] -= sl*omegai;
				}

				forAll(R.rhs(), s)
				{
				    label si = R.rhs()[s].index;
				    scalar sr = R.rhs()[s].stoichCoeff;
				    dcdt[si] += sr*omegai;
				}
			}

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
            }
        }
    }
}


template<class CompType, class ThermoType>
void Foam::ODEChemistryModel<CompType, ThermoType>::calculateElectronEnergySource()
{
    // The rate cache and the cacheable flags are set up by calculateWcf
    if
    (
        !this->chemistry_
     || cacheableRate_.size() != nReaction_
     || rateCacheCells_ != this->mesh().nCells()
    )
    {
        calculateWcf();
        return;
    }

    const volScalarField rho(this->thermo().rho());

    const labelList& kinds = reactionKinds();

    const scalarField& pCells = this->thermo().p();
    const scalarField& TCells = this->thermo().T();
    const scalarField& TeCells = this->thermo().Te();

    scalarField c(nSpecie_);

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    // The electron-impact reactions of calculateWcf, term by term in the
    // same order, so that the electron energy source is the same
    forAll(rho, celli)
    {
        eChemSource_[celli] = 0.0;
        dEChemSourceDTe_[celli] = 0.0;

        forAll(reactions_, i)
        {
            if (kinds[i] == electronElastic)
            {
                collFreq_[reactions_[i].lhs()[0].index][celli] = 0.0;
            }
        }

        const scalar rhoi = rho[celli];
        const scalar Ti = TCells[celli];
        const scalar pi = pCells[celli];
        const scalar Tei = TeCells[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
        }

        forAll(reactions_, i)
        {
            const label kind = kinds[i];

            if (kind != electronInelastic && kind != electronElastic)
            {
                continue;
            }

            const Reaction<ThermoType>& R = reactions_[i];

            scalar kfr, krr, kfp, krp, kfm, krm;
            cachedRate(i, celli, 0, Tei, pi, c, kfr, krr);
            cachedRate(i, celli, 1, Tei + 1e-08, pi, c, kfp, krp);
            cachedRate(i, celli, 2, Tei - 1e-08, pi, c, kfm, krm);

            const scalar omegai =
                krr == 0.0
              ? omegaForward(R, c, kfr)
              : omegaFromRates(R, c, kfr, krr, pf, cf, lRef, pr, cr, rRef);

            const scalar omegaiplus =
                krp == 0.0
              ? omegaForward(R, c, kfp)
              : omegaFromRates(R, c, kfp, krp, pf, cf, lRef, pr, cr, rRef);

            const scalar omegaiminus =
                krm == 0.0
              ? omegaForward(R, c, kfm)
              : omegaFromRates(R, c, kfm, krm, pf, cf, lRef, pr, cr, rRef);

            const scalar dOmegaiDTe = (omegaiplus - omegaiminus) / 2e-08;

            if (kind == electronInelastic)
            {
                scalar sl = R.lhs()[0].stoichCoeff;
                eChemSource_[celli] += sl*omegai*R.deltaE()*5.48579E-4;
                dEChemSourceDTe_[celli] += sl*dOmegaiDTe*R.deltaE()*5.48579E-4;
            }
            else
            {
                const label si = R.lhs()[1].index;
                const label sie = R.lhs()[0].index;
                eChemSource_[celli] += 3.0*omegai*(Tei-Ti)*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W();
                collFreq_[sie][celli] += omegai*5.48579E-4/Y_[sie][celli]/rhoi;
                dEChemSourceDTe_[celli] += 3.0*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W()*(dOmegaiDTe*(Tei-Ti) + omegai);
            }
        }
    }
}


template<class CompType, class ThermoType>
Foam::scalar Foam::ODEChemistryModel<CompType, ThermoType>::solve
(
    const scalar t0,
    const scalar deltaT
)
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    for (label i = 0; i < nSpecie_; i++)
    {
        RR_[i].setSize(rho.size());
    }

    if (!this->chemistry_)
    {
        return GREAT;
    }

    scalar deltaTMin = GREAT;

    tmp<volScalarField> thc = this->thermo().hc();
    const scalarField& hc = thc();

    forAll(rho, celli)
    {
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = 0.0;

        }

        scalar rhoi = rho[celli];
        scalar Ti = this->thermo().T()[celli];
        scalar hi = this->thermo().hs()[celli] + hc[celli];
        scalar pi = this->thermo().p()[celli];
        scalar Tei = this->thermo().Te()[celli];
        scalar Tioni = this->thermo().Tion()[celli];

        scalarField c(nSpecie_);
        scalarField c0(nSpecie_);
        scalarField dc(nSpecie_, 0.0);

        for (label i = 0; i < nSpecie_; i++)
        {
            c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
        }
        c0 = c;

        scalar t = t0;
        scalar tauC = this->deltaTChem_[celli];
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        // calculate the chemical source terms
        scalar cTot = 0.0;

        while (timeLeft > SMALL)
        {
            tauC = solver().solve(c, Ti, pi, Tei, Tioni, t, dt);
            t += dt;

            // update the temperature
            cTot = sum(c);
            ThermoType mixture(0.0*specieThermo_[0]);
            for (label i=0; i<nSpecie_; i++)
            {
                mixture += (c[i]/cTot)*specieThermo_[i];
            }
            Ti = mixture.TH(hi, Ti);

            timeLeft -= dt;
            this->deltaTChem_[celli] = tauC;
            dt = min(timeLeft, tauC);
            dt = max(dt, SMALL);
        }
        deltaTMin = min(tauC, deltaTMin);

        dc = c - c0;
        scalar WTot = 0.0;
        for (label i=0; i<nSpecie_; i++)
        {
            WTot += c[i]*specieThermo_[i].W();
        }
        WTot /= cTot;

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dc[i]*specieThermo_[i].W()/deltaT;
        }
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}


template<class CompType, class ThermoType>
inline Foam::scalarField&
Foam::ODEChemistryModel<CompType, ThermoType>::coeffs()
{
    return coeffs_;
}


template<class CompType, class ThermoType>
inline const Foam::scalarField&
Foam::ODEChemistryModel<CompType, ThermoType>::coeffs() const
{
    return coeffs_;
}


// ************************************************************************* //
