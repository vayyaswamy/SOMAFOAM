/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Class
Foam::electronThermalVelocity
Based on Thermal boundary conditions

SourceFiles
electronThermalVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef electronThermalVelocity_H
#define electronThermalVelocity_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class electronThermalVelocity
:
public fixedValueFvPatchVectorField
{
	// Private data
	scalar seec_;

	bool Edepend_;

    bool FE_; 

    scalar beta_;

    scalar wf_;
	
	public:

	//- Runtime type information
	TypeName("electronThermalVelocity");


	// Constructors

	//- Construct from patch and internal field
	electronThermalVelocity
	(
	 const fvPatch&,
	 const DimensionedField<vector, volMesh>&
	 );

	//- Construct from patch, internal field and dictionary
	electronThermalVelocity
	(
	 const fvPatch&,
	 const DimensionedField<vector, volMesh>&,
	 const dictionary&
	 );

	//- Construct by mapping given
	//  electronThermalVelocity onto a new patch
	electronThermalVelocity
	(
		const electronThermalVelocity&,
		const fvPatch&,
		const DimensionedField<vector, volMesh>&,
		const fvPatchFieldMapper&
	);

	//- Construct as copy
	electronThermalVelocity
	(
		const electronThermalVelocity&
	);

	//- Construct and return a clone
	virtual tmp<fvPatchVectorField> clone() const
	{
	  return tmp<fvPatchVectorField>
	(
	 new electronThermalVelocity(*this)
	 );
	}

	//- Construct as copy setting internal field reference
	electronThermalVelocity
	(
	 const electronThermalVelocity&,
	 const DimensionedField<vector, volMesh>&
	 );

	//- Construct and return a clone setting internal field reference
	virtual tmp<fvPatchVectorField> clone
	(
	 const DimensionedField<vector, volMesh>& iF
	) const
	{
		return tmp<fvPatchVectorField>
		(
			new electronThermalVelocity(*this, iF)
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
		const fvPatchVectorField&,
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
