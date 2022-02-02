/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::driftDiffusion

Description
    driftDiffusion equation - includes temperature in flux, fully implicit.
  
SourceFiles
    driftDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef driftDiffusion_H
#define driftDiffusion_H
#include "fvm.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class driftDiffusion Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class driftDiffusion
:
    public multiSpeciesPlasmaModel
{  

protected:
  
    // Protected data
  
        const PtrList<ThermoType>& speciesThermo_;

		volScalarField EON;

    // Protected member functions
    
        virtual inline void updateFlux
		(
			const label i,
			const volVectorField& E
		);

        virtual inline void transportCoeffInterpolate
		(
			const label i
		);

public:

		//- Runtime type information
		TypeName("driftDiffusion");
    
    // Constructors

        //- Construct from components
        driftDiffusion
        (
            hsCombustionThermo& thermo
        );

        //- Destructor
        ~driftDiffusion()
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
#   include "driftDiffusion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
