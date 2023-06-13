/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.


Class
Foam::ionHGLFlux
Based on boundary conditions described in Hagelaar et al.(2000) - https://journals.aps.org/pre/pdf/10.1103/PhysRevE.62.1452

SourceFiles
ionHGLFlux.C

\*---------------------------------------------------------------------------*/

#ifndef ionHGLFlux_H
#define ionHGLFlux_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ionHGLFlux
:
public fixedValueFvPatchVectorField
{
	// Private data
    word mobName_;
	scalar scalar1Data_;
	scalar scalar2Data_;
	scalar scalar3Data_;


	public:

	//- Runtime type information
	TypeName("ionHGLFlux");


	// Constructors

	//- Construct from patch and internal field
	ionHGLFlux
	(
	 const fvPatch&,
	 const DimensionedField<vector, volMesh>&
	 );

	//- Construct from patch, internal field and dictionary
	ionHGLFlux
	(
	 const fvPatch&,
	 const DimensionedField<vector, volMesh>&,
	 const dictionary&
	 );

	//- Construct by mapping given
	//  ionHGLFlux onto a new patch
	ionHGLFlux
	(
		const ionHGLFlux&,
		const fvPatch&,
		const DimensionedField<vector, volMesh>&,
		const fvPatchFieldMapper&
	);

	//- Construct as copy
	ionHGLFlux
	(
		const ionHGLFlux&
	);

	//- Construct and return a clone
	virtual tmp<fvPatchVectorField> clone() const
	{
	  return tmp<fvPatchVectorField>
	(
	 new ionHGLFlux(*this)
	 );
	}

	//- Construct as copy setting internal field reference
	ionHGLFlux
	(
	 const ionHGLFlux&,
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
			new ionHGLFlux(*this, iF)
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
