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
    dRRDi_(nSpecie_)
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
    Info<< "ODEChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
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
    scalarField c2(nSpecie_, 0.0);
    for (label i=0; i<nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    scalar kf, kr;

    //Info << "inside omega " << endl;

    if ((R.type() == "electronImpactInelasticArrheniusReaction") || (R.type() == "electronImpactInelasticTabularReaction"))
	{
      kf = R.kf(Te, p, c);
      kr = R.kr(kf, Te, p, c);
    }
    else if ((R.type() == "ionIrreversibleArrheniusReaction"))
	{
      kf = R.kf(Tion, p, c);
      kr = R.kr(kf, Tion, p, c);
    }
    else
	{
      kf = R.kf(T, p, c);
      kr = R.kr(kf, T, p, c);
    }

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
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);


    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > ROOTVSMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
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
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        scalar exp = R.rhs()[srRef].exponent;
        if (exp<1.0)
        {
            if (cr>ROOTVSMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
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
                            kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c2[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c2[si], el);
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
                            kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c2[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= pow(c2[si], er);
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
    }

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
                                    kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                                }
                                else
                                {
                                    kf = 0.0;
                                }
                            }
                            else
                            {
                                kf *= el*pow(c2[si], el - 1.0);
                            }
                        }
                        else
                        {
                            kf *= pow(c2[si], el);
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
                                    kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                                }
                                else
                                {
                                    kr = 0.0;
                                }
                            }
                            else
                            {
                                kr *= er*pow(c2[si], er - 1.0);
                            }
                        }
                        else
                        {
                            kr *= pow(c2[si], er);
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
    }

    if (this->chemistry_)
    {
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
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];
            scalar Tei = this->thermo().Te()[celli];
            scalar Tioni = this->thermo().Tion()[celli];

            scalarField c(nSpecie_);
            scalarField dcdt(nEqns(), 0.0);
            scalarField c2(nSpecie_, 0.0);

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
                        if (si == sj)
                        {
                            if (el < 1.0)
                            {
                                if (c2[si]>ROOTVSMALL)
                                {
                                    kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                                }
                                else
                                {
                                    kf = 0.0;
                                }
                            }
                            else
                            {
                                kf *= el*pow(c2[si], el - 1.0);
                            }
                        }
                        else
                        {
                            kf *= pow(c2[si], el);
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
                                    kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                                }
                                else
                                {
                                    kr = 0.0;
                                }
                            }
                            else
                            {
                                kr *= er*pow(c2[si], er - 1.0);
                            }
                        }
                        else
                        {
                            kr *= pow(c2[si], er);
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
					const label sie = R.lhs()[0].index; 
					eChemSource_[celli] += 3.0*omegai*(Tei-Ti)*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W(); //electron weight square
					collFreq_[sie][celli] += omegai*5.48579E-4/Y_[sie][celli]/rhoi;
                    dEChemSourceDTe_[celli] += 3.0*3.00938919241e-7*1.38064852E-23/specieThermo_[si].W()*(dOmegaiDTe*(Tei-Ti) + omegai);
				  // need to multiply by electron mass and divide by other species mass
				}
				if ((R.type() == "ionElasticArrheniusReaction"))
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
