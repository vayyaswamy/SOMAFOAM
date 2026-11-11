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

    driveTable_()
{
    if (found("zeroDDrive"))
    {
        const dictionary& driveDict = subDict("zeroDDrive");

        driveMode_ = word(driveDict.lookup("mode"));

        if
        (
            driveMode_ != "none"
         && driveMode_ != "reducedField"
         && driveMode_ != "electronTemperature"
        )
        {
            FatalIOErrorIn("zeroD::zeroD", driveDict)
                << "unknown zeroDDrive mode " << driveMode_
                << ", use none, reducedField or electronTemperature"
                << exit(FatalIOError);
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::zeroD<ThermoType>::input()
{
    // The species dictionaries (mobility models among them) are read here
    multiSpeciesPlasmaModel::input();

    if
    (
        driveMode_ == "reducedField"
     && mobilityModel_[eIndex_] != "eTemp"
     && mobilityModel_[eIndex_] != "constant"
    )
    {
        FatalErrorIn("zeroD::input")
            << "reducedField drive needs the electron mobilityModel"
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
