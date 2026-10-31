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

// Whether the abscissae of a table increase strictly, remembered for the
// last few tables (the transport tables are fixed for the whole run)
inline bool plasmaInterpolateIncreasing(const scalarField& xOld)
{
    static const label nMemo = 16;
    static const scalar* ptr[nMemo] = {};
    static label size[nMemo] = {};
    static bool increasing[nMemo] = {};
    static label next = 0;

    for (label k = 0; k < nMemo; ++k)
    {
        if (ptr[k] == xOld.cdata() && size[k] == xOld.size())
        {
            return increasing[k];
        }
    }

    bool inc = true;

    for (label i = 1; i < xOld.size(); ++i)
    {
        if (!(xOld[i] > xOld[i-1]))
        {
            inc = false;
            break;
        }
    }

    ptr[next] = xOld.cdata();
    size[next] = xOld.size();
    increasing[next] = inc;
    next = (next + 1) % nMemo;

    return inc;
}


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

    // Strictly increasing abscissae (the transport tables): bisection for
    // the interval the scans below find, with the same arithmetic. lo is
    // the last row with xOld[lo] <= x, hi the first with xOld[hi] >= x
    if (n > 1 && x == x && plasmaInterpolateIncreasing(xOld))
    {
        if (x < xOld[0])
        {
            return yOld[0]
                + ((x - xOld[0])/(xOld[1] - xOld[0]))*(yOld[1] - yOld[0]);
        }

        if (x >= xOld[n-1])
        {
            return yOld[n-1];
        }

        label l = 0;
        label h = n - 1;

        while (h - l > 1)
        {
            const label mid = (l + h)/2;

            if (xOld[mid] <= x)
            {
                l = mid;
            }
            else
            {
                h = mid;
            }
        }

        if (xOld[l] == x)
        {
            return yOld[l];
        }

        return yOld[l]
            + ((x - xOld[l])/(xOld[h] - xOld[l]))*(yOld[h] - yOld[l]);
    }

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
