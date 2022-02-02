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
externalCircuitPotential.C

\*---------------------------------------------------------------------------*/

#ifndef externalCircuitPotential_H
#define externalCircuitPotential_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class externalCircuitPotential
:
public fixedValueFvPatchScalarField
{
    // Private data
    scalarField amplitude_;
    word modelName_;
    scalar frequency_;
    scalar bias_;
    scalar R_;
    //scalar iSum_;
    //scalar iRMS_;

    //scalar iDesired_


    	
  public:

    //- Runtime type information
    TypeName("externalCircuitPotential");


    // Constructors

    //- Construct from patch and internal field
    externalCircuitPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    externalCircuitPotential
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  plasmaPotential onto a new patch
    externalCircuitPotential
    (
     const externalCircuitPotential&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    externalCircuitPotential
    (
     	const externalCircuitPotential&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
      	return tmp<fvPatchScalarField>
		(
		 new externalCircuitPotential(*this)
		);
    }

    //- Construct as copy setting internal field reference
    externalCircuitPotential
    (
		 const externalCircuitPotential&,
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
			new externalCircuitPotential(*this, iF)
		);
    }


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
