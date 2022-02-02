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

inline Foam::ModifiedArrheniusReactionRate::ModifiedArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar omega
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    omega_(omega)
{}


inline Foam::ModifiedArrheniusReactionRate::ModifiedArrheniusReactionRate
(
    const speciesTable&,
    Istream& is
)
:
    A_(readScalar(is.readBegin("ModifiedArrheniusReactionRate(Istream&)"))),
    beta_(readScalar(is)),
    Ta_(readScalar(is)),
    omega_(readScalar(is))
{
    is.readEnd("ModifiedArrheniusReactionRate(Istream&)");
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ModifiedArrheniusReactionRate::operator()
(
    const scalar T,
    const scalar,
    const scalarField&
) const
{
    scalar ak = A_;

    if (mag(beta_) > VSMALL)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > VSMALL)
    {
      scalar exparg = pow(Ta_/T,omega_);
      ak *= exp(-exparg);
    }

    return ak;
}

inline Ostream& operator<<
(
    Ostream& os,
    const ModifiedArrheniusReactionRate& arr
)
{
    os  << token::BEGIN_LIST
        << arr.A_ << token::SPACE << arr.beta_ << token::SPACE << arr.Ta_
        << token::SPACE << arr.omega_ << token::END_LIST;
    return os;
}

}
// ************************************************************************* //
