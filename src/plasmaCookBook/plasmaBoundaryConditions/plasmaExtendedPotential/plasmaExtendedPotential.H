/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
Foam::plasmaExtendedPotential

Description
    Abstract base class for generalized electric potential boundary conditions.

SourceFiles
plasmaExtendedPotential.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaExtendedPotential_H
#define plasmaExtendedPotential_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaExtendedPotential
:
public fixedValueFvPatchScalarField
{
    // Private data

    scalarField voltage() const;
	
  public:

    //- Runtime type information
    TypeName("plasmaExtendedPotential");


    // Constructors

    //- Construct from patch and internal field
    plasmaExtendedPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    plasmaExtendedPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  plasmaExtendedPotential onto a new patch
    plasmaExtendedPotential
    (
     const plasmaExtendedPotential&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    plasmaExtendedPotential
    (
     	const plasmaExtendedPotential&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
      	return tmp<fvPatchScalarField>
		(
		 new plasmaExtendedPotential(*this)
		);
    }

    //- Construct as copy setting internal field reference
    plasmaExtendedPotential
    (
		 const plasmaExtendedPotential&,
		 const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
     	const DimensionedField<scalar, volMesh>& iF
    ) const
    {
		return tmp<fvPatchScalarField>
		(
			new plasmaExtendedPotential(*this, iF)
		);
    }

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMap
    (
      const fvPatchFieldMapper&
    );

    //- Reverse map the given fvPatchField onto this fvPatchField
    virtual void rmap
    (
      const fvPatchScalarField&,
      const labelList&
    );
    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
