/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "electronVelocity.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

//#include "fvPatchFieldMapper.H"
//#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electronVelocity::
electronVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Te_(0.0)
{}


Foam::electronVelocity::
electronVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    Te_(readScalar(dict.lookup("Te")))
{}


Foam::electronVelocity::
electronVelocity
(
    const electronVelocity& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Te_(ptf.Te_)
{}


Foam::electronVelocity::
electronVelocity
(
    const electronVelocity& eFpvf
)
:
    fixedValueFvPatchVectorField(eFpvf),
    Te_(eFpvf.Te_)
{}


Foam::electronVelocity::
electronVelocity
(
    const electronVelocity& eFpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(eFpvf, iF),
    Te_(eFpvf.Te_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electronVelocity::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // const fvMesh& mesh = internalField().mesh();
    
    const dictionary& physicalProperties = db().lookupObject<IOdictionary>
	    (
	        "physicalProperties"
	    );
    dimensionedScalar me(physicalProperties.lookup("me"));

    scalar e = 1.602e-19;

    scalar pi = acos(-1.0);

//    const fvPatchField<scalar>& nef =
//	    patch().lookupPatchField<volScalarField, scalar>("ne");

    //const fvPatchField<vector>& Ef =
//	    patch().lookupPatchField<volVectorField, vector>("E");
  
  //  const volVectorField& mome = 
//	    db().objectRegistry::lookupObject<volVectorField>("mome");
   
    vectorField n = patch().nf();
    
    vectorField Uef = 0.0*n; 

    forAll(Uef, facei)
    {
 //       label faceCelli = patch().faceCells()[facei];

	Uef[facei] = 0.25*sqrt(8*e*Te_/pi/me.value())*n[facei];
    }

    operator == (Uef);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::electronVelocity::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Te")
	    << Te_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        electronVelocity
    );
}

// ************************************************************************* //
