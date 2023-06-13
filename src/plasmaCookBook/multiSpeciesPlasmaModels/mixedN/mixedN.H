/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::mixedN

Description
    mixedN
  
SourceFiles
    mixedN.C

\*---------------------------------------------------------------------------*/

#ifndef mixedN_H
#define mixedN_H
#include "fvm.H"
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class mixed Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class mixedN
:
    public multiSpeciesPlasmaModel
{  

protected:
  
    // Protected data
  
        const PtrList<ThermoType>& speciesThermo_;

		scalar momentumCount;

		scalar ddCount;

		volScalarField EON;

    // Protected member functions

        virtual void mixedModelInput();

        virtual inline void updateVelocity
		(
			const psiChemistryModel& chemistry,
			const volVectorField& E,
			const label i
		);

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
		TypeName("mixedN");
    
    // Constructors

        //- Construct from components
        mixedN
        (
            hsCombustionThermo& thermo
        );

        //- Destructor
        ~mixedN()
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
#   include "mixedN.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
