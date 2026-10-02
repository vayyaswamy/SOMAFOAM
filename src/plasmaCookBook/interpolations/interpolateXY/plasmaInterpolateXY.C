/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.
\*---------------------------------------------------------------------------*/

#include "plasmaInterpolateXY.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Field<Type> plasmaInterpolateXY
(
    const scalarField& xNew,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    scalarField yNew(xNew.size());


    forAll(xNew, i)
    {
        yNew[i] = plasmaInterpolateXY(xNew[i], xOld, yOld);
    }

    return yNew;
}


template<class Type>
Type plasmaInterpolateXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    label n = xOld.size();

    label lo = 0;
    for (lo=0; lo<n && xOld[lo]>x; ++lo)
    {}

    label low = lo;
    if (low < n)
    {
        for (label i=low; i<n; ++i)
        {
            if (xOld[i] > xOld[lo] && xOld[i] <= x)
            {
                lo = i;
            }
        }
    }

    label hi = 0;
    for (hi=0; hi<n && xOld[hi]<x; ++hi)
    {}

    label high = hi;
    if (high < n)
    {
        for (label i=high; i<n; ++i)
        {
            if (xOld[i] < xOld[hi] && xOld[i] >= x)
            {
                hi = i;
            }
        }
    }


    if (lo<n && hi<n && lo != hi)
    {
        return yOld[lo]
            + ((x - xOld[lo])/(xOld[hi] - xOld[lo]))*(yOld[hi] - yOld[lo]);
    }
    else if (lo == hi)
    {
        return yOld[lo];
    }
    else if (lo == n)
    {
        return yOld[hi]
            + ((x - xOld[hi])/(xOld[hi+1] - xOld[hi]))*(yOld[hi+1] - yOld[hi]);
    }
    else
    {
        return yOld[lo];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
