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

    //Info << "xNew = " << xNew << endl;
    //Info << "xOld = " << xOld << endl;
    //Info << "yOld = " << yOld << endl;

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
        //return yOld[lo];
        //Info << "lo = " << lo << endl;
        //Info << "n = " << n << endl;
        //Info << "value = " << x << " " << xOld[lo] << " " << xOld[lo+1] << endl;
        //Info << "value = " << yOld[lo] << " " << yOld[lo+1] << endl;
        //Info << exp(log(yOld[lo]) + 
        //    (log(x) - log(xOld[lo]))/(log(xOld[lo+1]) - log(xOld[lo]))*(log(yOld[lo+1]) - log(yOld[lo]))); 
        //return exp(log(yOld[lo]) + 
        //    (log(x) - log(xOld[lo]))/(log(xOld[lo+1]) - log(xOld[lo]))*(log(yOld[lo+1]) - log(yOld[lo])));
        return yOld[lo];
    }
    else if (lo == n)
    {
        return yOld[hi]
            + ((x - xOld[hi])/(xOld[hi+1] - xOld[hi]))*(yOld[hi+1] - yOld[hi]);
    }
    else
    {
        //Info << "x = " << x << endl;
        //Info << "lo = " << lo << endl;
        //Info << "hi = " << hi << endl;
        //Info <<"n = " << n << endl;
        //return exp(log(yOld[lo]) + 
        //    (log(x) - log(xOld[lo]))/(log(xOld[lo+1]) - log(xOld[lo]))*(log(yOld[lo+1]) - log(yOld[lo])));
        return yOld[lo];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
