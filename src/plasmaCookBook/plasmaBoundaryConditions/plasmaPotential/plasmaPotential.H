/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
Foam::plasmaPotential

Description
    Abstract base class for electric potential boundary conditions.

SourceFiles
plasmaPotential.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaPotential_H
#define plasmaPotential_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaPotential
:
public fixedValueFvPatchScalarField
{
    // Private data
    scalarField amplitude_;
    word modelName_;
    scalar frequency_;
    scalar bias_;
    	
  public:

    //- Runtime type information
    TypeName("plasmaPotential");


    // Constructors

    //- Construct from patch and internal field
    plasmaPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    plasmaPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  plasmaPotential onto a new patch
    plasmaPotential
    (
     const plasmaPotential&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    plasmaPotential
    (
     	const plasmaPotential&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
      	return tmp<fvPatchScalarField>
		(
		 new plasmaPotential(*this)
		);
    }

    //- Construct as copy setting internal field reference
    plasmaPotential
    (
		 const plasmaPotential&,
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
			new plasmaPotential(*this, iF)
		);
    }


    // Member functions

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
