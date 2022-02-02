/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Class
Foam::electronThermalFlux
Based on Thermal boundary conditions

SourceFiles
electronThermalFlux.C

\*---------------------------------------------------------------------------*/

#ifndef electronThermalFlux_H
#define electronThermalFlux_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class electronThermalFlux
:
public fixedValueFvPatchVectorField
{
	// Private data
	scalar scalar1Data_;
	
	public:

	//- Runtime type information
	TypeName("electronThermalFlux");


	// Constructors

	//- Construct from patch and internal field
	electronThermalFlux
	(
	 const fvPatch&,
	 const DimensionedField<vector, volMesh>&
	 );

	//- Construct from patch, internal field and dictionary
	electronThermalFlux
	(
	 const fvPatch&,
	 const DimensionedField<vector, volMesh>&,
	 const dictionary&
	 );

	//- Construct by mapping given
	//  electronThermalFlux onto a new patch
	electronThermalFlux
	(
		const electronThermalFlux&,
		const fvPatch&,
		const DimensionedField<vector, volMesh>&,
		const fvPatchFieldMapper&
	);

	//- Construct as copy
	electronThermalFlux
	(
		const electronThermalFlux&
	);

	//- Construct and return a clone
	virtual tmp<fvPatchVectorField> clone() const
	{
	  return tmp<fvPatchVectorField>
	(
	 new electronThermalFlux(*this)
	 );
	}

	//- Construct as copy setting internal field reference
	electronThermalFlux
	(
	 const electronThermalFlux&,
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
			new electronThermalFlux(*this, iF)
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
