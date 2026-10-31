/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors
License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.
(This code was written by Abhishek Kumar Verma a former PhD student at UC Merced)
Class
Foam::plasmaPower
Description
    Abstract base class for electric potential boundary conditions.
SourceFiles
plasmaPower.C
\*---------------------------------------------------------------------------*/

#ifndef plasmaPower_H
#define plasmaPower_H
//#include <iostream>
#include "fixedValueFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaPower
:
public fixedValueFvPatchScalarField
{
    // Private data
    scalar currentTimee=this->db().time().value() ;
    scalar previousTimee=this->db().time().value();
    //scalarField amplitude_;
    scalar amplitude_;
    //word modelName_;
    scalar frequency_;
    scalar nCycles_;
    scalar power_;
    scalar bias_;
    scalar sumPowerr_=0;
    int nCountter=0;
    int timesCorrect=0;
    scalar averagePower = 0;
    scalar dampingFactor_;

    OFstream* powerLogFilePtr_;
    	
  public:

    //- Runtime type information
    TypeName("plasmaPower");


    // Constructors

    //- Construct from patch and internal field
    plasmaPower
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&
     );

    //- Construct from patch, internal field and dictionary
    plasmaPower
    (
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const dictionary&
     );

    //- Construct by mapping given
    //  plasmaPower onto a new patch
    plasmaPower
    (
     const plasmaPower&,
     const fvPatch&,
     const DimensionedField<scalar, volMesh>&,
     const fvPatchFieldMapper&
     );

    //- Construct as copy
    plasmaPower
    (
     	const plasmaPower&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
      	return tmp<fvPatchScalarField>
		(
		 new plasmaPower(*this)
		);
    }

    //- Construct as copy setting internal field reference
    plasmaPower
    (
		 const plasmaPower&,
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
			new plasmaPower(*this, iF)
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