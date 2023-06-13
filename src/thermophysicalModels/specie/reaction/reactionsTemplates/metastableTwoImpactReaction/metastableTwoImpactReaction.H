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
    Foam::metastableTwoImpactReaction

Description
    Simple extension of Reaction to handle electron impact inelastic reactions

SourceFiles
    metastableTwoImpactReaction.C

\*---------------------------------------------------------------------------*/

#ifndef metastableTwoImpactReaction_H
#define metastableTwoImpactReaction_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class metastableTwoImpactReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ReactionRate>
class metastableTwoImpactReaction
:
    public Reaction<ReactionThermo>
{
    // Private data
    ReactionRate k_;
    scalar deltaE_;
    scalar fractional_;

    // Private Member Functions

    //- Disallow default bitwise assignment
    void operator=
    (
        const metastableTwoImpactReaction<ReactionThermo, ReactionRate>&
    );


public:

    //- Runtime type information
    TypeName("metastableTwoImpact");


    // Constructors

    //- Construct from components
    metastableTwoImpactReaction
    (
        const Reaction<ReactionThermo>& reaction,
        const ReactionRate& reactionRate,
    	const scalar deltaE, 
        const scalar fractional
    );

    //- Construct as copy given new speciesTable
    metastableTwoImpactReaction
    (
        const metastableTwoImpactReaction<ReactionThermo, ReactionRate>&,
        const speciesTable& species
    );

    //- Construct from Istream
    metastableTwoImpactReaction
    (
        const speciesTable& species,
        const HashPtrTable<ReactionThermo>& thermoDatabase,
        Istream& is
    );

    //- Construct and return a clone
    virtual autoPtr<Reaction<ReactionThermo> > clone() const
    {
        return autoPtr<Reaction<ReactionThermo> >
        (
            new metastableTwoImpactReaction<ReactionThermo, ReactionRate>(*this)
        );
    }

    //- Construct and return a clone with new speciesTable
    virtual autoPtr<Reaction<ReactionThermo> > clone
    (
        const speciesTable& species
    ) const
    {
        return autoPtr<Reaction<ReactionThermo> >
        (
            new metastableTwoImpactReaction<ReactionThermo, ReactionRate>
            (
                *this,
                species
            )
        );
    }



    // Member Functions

    //- Destructor
    virtual ~metastableTwoImpactReaction()
    {

    }
    
    // metastableTwoImpactReaction rate coefficients

    //- Forward rate constant
    virtual scalar kf
    (
        const scalar T,
        const scalar p,
        const scalarField& c
    ) const;

    //- deltaE
    virtual scalar deltaE() const;

    //- fractional
    virtual scalar fractional() const; 

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "metastableTwoImpactReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
