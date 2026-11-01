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

Class
    Foam::plasmaFieldChemistryModel

Description
    Plasma chemistry evaluated reaction by reaction over all cells, for the
    direct solution of the plasma equations (plasmaChemistryModel
    directSolution). It reads the same reactions as ODEChemistryModel and
    gives the same source terms, rates, linearization, electron energy
    source and collision frequencies (calculateWcf), with three
    differences in how they are found:

    - each rate coefficient is evaluated once per cell and call, for the
      whole field, instead of once for the rate and twice more for its
      electron temperature derivative;
    - tabulated electron-impact rates are interpolated directly (ln k or,
      with rateInterpolation linear, k linearly in Te), with the derivative
      of the table segment instead of a central difference;
    - the concentrations, stoichiometry and kinds of the reactions are
      held in plain arrays, without the reaction objects in the cell loops;
    - the rate coefficients of a reaction are kept while the temperature
      and pressure they depend on do not change (the species updates of a
      time step leave Te unchanged).

    Reactions with reverse rates are not handled this way: a mechanism
    that has any falls back to ODEChemistryModel::calculateWcf.

    Selected in constant/chemistryProperties with
    \verbatim
    psiChemistryModel  plasmaFieldChemistryModel<constGasThermoPhysics>;
    \endverbatim

SourceFiles
    plasmaFieldChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaFieldChemistryModel_H
#define plasmaFieldChemistryModel_H

#include "ODEChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class plasmaFieldChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class plasmaFieldChemistryModel
:
    public ODEChemistryModel<CompType, ThermoType>
{
    // Private data

        //- Set up from the reactions on the first call
        bool setUp_;

        //- Some reaction has a reverse rate: use ODEChemistryModel
        bool reversible_;

        //- Per reaction: species, stoichiometric coefficients and exponents
        //  of the left- and right-hand sides
        List<labelList> lhsSpecie_;
        List<scalarList> lhsCoeff_;
        List<scalarList> lhsExp_;
        List<labelList> rhsSpecie_;
        List<scalarList> rhsCoeff_;

        //- Per reaction and left-hand side species: net production
        //  coefficient of that species in the reaction (rhs - lhs)
        List<scalarList> lhsNet_;

        //- Per reaction: ln k against Te of a tabulated electron-impact
        //  rate with increasing abscissae (empty otherwise)
        List<scalarList> tableT_;
        List<scalarList> tableLnK_;

        //- Per reaction: the rate does not depend on the concentrations
        boolList fieldRate_;

        //- Per reaction (rates that do not depend on the concentrations):
        //  the temperatures and pressures at which kf and dkf/dTe were last
        //  evaluated, and their values. The species updates of a time step
        //  do not change Te, and their rates are taken from here
        List<scalarField> cacheT_;
        List<scalarField> cacheP_;
        List<scalarField> cacheKf_;
        List<scalarField> cacheDkf_;
        boolList cacheDerivative_;

        //- Work fields over the cells
        scalarField kf_;
        scalarField dkfdTe_;
        scalarField prodC_;
        PtrList<scalarField> c_;
        PtrList<scalarField> dcdt_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        plasmaFieldChemistryModel(const plasmaFieldChemistryModel&);
        void operator=(const plasmaFieldChemistryModel&);

        //- Read the stoichiometry and the rate tables of the reactions
        void setUp();

        //- Concentrations of all species [kmol/m3], and work fields
        void concentrations(const scalarField& rho);

        //- Forward rate coefficient of reaction i in every cell at the
        //  temperatures Tr, and (withDerivative) its derivative in Tr
        void forwardRate
        (
            const label i,
            const scalarField& Tr,
            const scalarField& p,
            const bool withDerivative
        );

        //- Product of the concentrations of the left-hand side of reaction i
        //  (prodC_), omega = kf prodC
        void reactantProduct(const label i);

        //- Electron energy source and collision frequency of the
        //  electron-impact reaction i, kf_, dkfdTe_ and prodC_ set
        void addElectronTerms
        (
            const label i,
            const label kind,
            const scalarField& rho,
            const scalarField& T,
            const scalarField& Te
        );


public:

    //- Runtime type information
    TypeName("plasmaFieldChemistryModel");


    // Constructors

        //- Construct from components
        plasmaFieldChemistryModel
        (
            const fvMesh& mesh,
            const objectRegistry& obj,
            const word& compTypeName,
            const word& thermoTypeName
        );


    //- Destructor
    virtual ~plasmaFieldChemistryModel();


    // Member Functions

        //- Source terms, linearization, electron energy source and
        //  collision frequencies
        virtual void calculateWcf();

        //- Only the electron energy source, its Te derivative and the
        //  electron collision frequency
        virtual void calculateElectronEnergySource();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "plasmaFieldChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
