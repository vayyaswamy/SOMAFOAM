/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::momentum

Description
    momentum equation.

    Each charged specie has its own momentum equation; the continuity
    equations are solved for the number densities N (the quantity held in
    the 0/<specie> files and assumed by the plasma boundary conditions).
    Neutral species use drift-free diffusion.

SourceFiles
    momentum.C

\*---------------------------------------------------------------------------*/

#ifndef momentum_H
#define momentum_H
#include "fvm.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class momentum Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class momentum
:
    public multiSpeciesPlasmaModel
{

protected:

        const PtrList<ThermoType>& speciesThermo_;

        volScalarField EON;

    // Protected member functions

        virtual inline void updateVelocity
        (
            const psiChemistryModel& chemistry,
            const volVectorField& E,
            const label i
        );

        virtual inline void transportCoeffInterpolate
        (
            const label i
        );

        //- Solve the number density equation of specie i
        void solveSpecie
        (
            const label i,
            psiChemistryModel& chemistry,
            const volVectorField& E
        );

public:

    //- Runtime type information
    TypeName("momentum");

    // Constructors

        //- Construct from components
        momentum
        (
            hsCombustionThermo& thermo
        );

      //- Destructor
      ~momentum()
      {}

    // Member Functions

        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }

        virtual inline scalar correct
        (
            psiChemistryModel& chemistry,
            const volVectorField& E,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );

        virtual bool read();
  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "momentum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
