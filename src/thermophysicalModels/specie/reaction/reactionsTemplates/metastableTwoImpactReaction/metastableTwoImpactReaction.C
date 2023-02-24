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

#include "metastableTwoImpactReaction.H"

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
metastableTwoImpactReaction<ReactionThermo, ReactionRate>::metastableTwoImpactReaction
(
    const Reaction<ReactionThermo>& reaction,
    const ReactionRate& k,
    const scalar deltaE,
    const scalar fractional
)
:
    Reaction<ReactionThermo>(reaction),
    k_(k),
    deltaE_(deltaE),
    fractional_(fractional)
{}


template<class ReactionThermo, class ReactionRate>
metastableTwoImpactReaction<ReactionThermo, ReactionRate>::
metastableTwoImpactReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, is),
    k_(species, is),
    deltaE_(readScalar(is.readBegin("metastableTwoImpactReaction(Istream&)")))
{
    is.readEnd("metastableTwoImpactReaction(Istream&)");

    fractional_=readScalar(is.readBegin("metastableTwoImpactReaction(Istream&)"));
    is.readEnd("metastableTwoImpactReaction(Istream&)");
}


template<class ReactionThermo, class ReactionRate>
metastableTwoImpactReaction<ReactionThermo, ReactionRate>::
metastableTwoImpactReaction
(
    const metastableTwoImpactReaction<ReactionThermo,ReactionRate>& irr,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>(irr, species),
    k_(irr.k_),
    deltaE_(irr.deltaE_),
    fractional_(irr.fractional_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
scalar metastableTwoImpactReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return k_(T, p, c);
}

template<class ReactionThermo, class ReactionRate>
scalar metastableTwoImpactReaction<ReactionThermo, ReactionRate>::deltaE() const
{
    return deltaE_;
}

template<class ReactionThermo, class ReactionRate>
scalar metastableTwoImpactReaction<ReactionThermo, ReactionRate>::fractional() const
{
    return fractional_;
}

template<class ReactionThermo, class ReactionRate>
void metastableTwoImpactReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);
    os  << token::SPACE << k_;
    os.writeKeyword("deltaE") << deltaE_ << token::END_STATEMENT << nl;
    os.writeKeyword("fractional") << fractional_ << token::END_STATEMENT << nl;
}

}
// ************************************************************************* //
