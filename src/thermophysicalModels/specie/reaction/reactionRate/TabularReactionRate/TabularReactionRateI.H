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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from components
inline TabularReactionRate::TabularReactionRate
(
 const fileName& filename
)
:
  filename_(filename),
  lookup_(filename_),
  increasing_(increasingRows(lookup_))
{}

//- Construct from Istream
inline TabularReactionRate::TabularReactionRate
(
    const speciesTable&,
    Istream& is
)
:
  filename_((is.readBegin("TabularReactionRate(Istream&)"))),
  lookup_(filename_),
  increasing_(increasingRows(lookup_))
{
    is.readEnd("TabularReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar TabularReactionRate::operator()
(
    const scalar T,
    const scalar,
    const scalarField&
) const
{
    scalar ak;

    const label n = lookup_.size();

    if (increasing_ && n > 1 && T >= lookup_[0].first() && T < lookup_[n-1].first())
    {
        // lookup_[lo].first() <= T < lookup_[hi].first(), hi = lo + 1: the
        // interval of interpolationTable::operator(), which searches it
        // linearly from the first row
        label lo = 0;
        label hi = n - 1;

        while (hi - lo > 1)
        {
            const label mid = (lo + hi)/2;

            if (T >= lookup_[mid].first())
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }

        // The interpolation of interpolationTable::operator(), term by term
        ak =
            lookup_[lo].second()
          + (lookup_[hi].second() - lookup_[lo].second())
           *(T - lookup_[lo].first())
           /(lookup_[hi].first() - lookup_[lo].first());
    }
    else
    {
        ak = lookup_(T);
    }

    return exp(ak);
}

inline Ostream& operator<<
(
    Ostream& os,
    const TabularReactionRate& arr
)
{
    os  << token::BEGIN_LIST
        << arr.filename_ << token::END_LIST;
    return os;
}

}
// ************************************************************************* //
