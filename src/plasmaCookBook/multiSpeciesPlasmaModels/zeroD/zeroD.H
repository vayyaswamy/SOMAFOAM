/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Class
    Foam::zeroD

Description
    zeroD equation.
  
SourceFiles
    zeroD.C

\*---------------------------------------------------------------------------*/

#ifndef zeroD_H
#define zeroD_H
#include "fvm.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class zeroD Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class zeroD
:
    public multiSpeciesPlasmaModel
{  

protected:
  
    // Protected data
  
        const PtrList<ThermoType>& speciesThermo_;

public:

		//- Runtime type information
		TypeName("zeroD");
    
    // Constructors

        //- Construct from components
        zeroD
        (
            hsCombustionThermo& thermo
        );

        //- Destructor
        ~zeroD()
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
                
        //- Read dictionary
        virtual bool read();
  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "zeroD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
