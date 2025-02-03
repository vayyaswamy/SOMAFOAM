/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
Foam::currentControlledElectrode

Description
    Abstract base class for electric potential boundary conditions.

SourceFiles
currentControlledElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef currentControlledElectrode_H
#define currentControlledElectrode_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class currentControlledElectrode
:
public fixedGradientFvPatchScalarField
{
    // Private data
    scalarField amplitude_;
    word modelName_;
    scalar frequency_;
    	
  public:

    //- Runtime type information
    TypeName("currentControlledElectrode");


    // Constructors

    //- Construct from patch and internal field
    currentControlledElectrode
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    currentControlledElectrode
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  currentControlledElectrode onto a new patch
    currentControlledElectrode
    (
     const currentControlledElectrode&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    currentControlledElectrode
    (
     	const currentControlledElectrode&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
      	return tmp<fvPatchScalarField>
		(
		 new currentControlledElectrode(*this)
		);
    }

    //- Construct as copy setting internal field reference
    currentControlledElectrode
    (
		 const currentControlledElectrode&,
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
			new currentControlledElectrode(*this, iF)
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
