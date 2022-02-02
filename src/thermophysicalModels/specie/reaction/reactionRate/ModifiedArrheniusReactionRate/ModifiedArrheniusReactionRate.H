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
    Foam::ModifiedArrheniusReactionRate

Description
    Modified Arrhenius reaction rate given by:

        k = A * T^beta * exp(-(Ta/T)^omega)

SourceFiles
    ModifiedArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef ModifiedArrheniusReactionRate_H
#define ModifiedArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ModifiedArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class ModifiedArrheniusReactionRate
{
    // Private data

        scalar A_;
        scalar beta_;
        scalar Ta_;
        scalar omega_;

public:

    // Constructors

        //- Construct from components
        inline ModifiedArrheniusReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const scalar omega
        );

        //- Construct from Istream
        inline ModifiedArrheniusReactionRate
        (
            const speciesTable& species,
            Istream& is
        );

    // Member Functions

        //- Return the type name
        static word type()
        {
            return "ModifiedArrhenius";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar p,
            const scalarField& c
        ) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const ModifiedArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ModifiedArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
