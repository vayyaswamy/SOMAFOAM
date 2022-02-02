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

#include "ionVelocity.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionVelocity::
ionVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ti_(0.0),
    mi_(0.0),
    variable_("none")
{}


Foam::ionVelocity::
ionVelocity
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    Ti_(readScalar(dict.lookup("Ti"))),
    mi_(readScalar(dict.lookup("mi"))),
    variable_(word(dict.lookup("variable")))
{}


Foam::ionVelocity::
ionVelocity
(
    const ionVelocity& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Ti_(ptf.Ti_),
    mi_(ptf.mi_),
    variable_(ptf.variable_)
{}


Foam::ionVelocity::
ionVelocity
(
    const ionVelocity& eFpvf
)
:
    fixedValueFvPatchVectorField(eFpvf),
    Ti_(eFpvf.Ti_),
    mi_(eFpvf.mi_),
    variable_(eFpvf.variable_)
{}


Foam::ionVelocity::
ionVelocity
(
    const ionVelocity& eFpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(eFpvf, iF),
    Ti_(eFpvf.Ti_),
    mi_(eFpvf.mi_),
    variable_(eFpvf.variable_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ionVelocity::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // const fvMesh& mesh = internalField().mesh();
    

    scalar e = 1.602e-19;

    scalar pi = acos(-1.0);

//    const fvPatchField<scalar>& nef =
//	    patch().lookupPatchField<volScalarField, scalar>("ne");

    //const fvPatchField<vector>& Ef =
//	    patch().lookupPatchField<volVectorField, vector>("E");
  
    const volVectorField& Ui = 
	    db().objectRegistry::lookupObject<volVectorField>(variable_);
   
    vectorField n = patch().nf();
    
    vectorField Uif = 0.0*n; 

    forAll(Uif, facei)
    {
        label faceCelli = patch().faceCells()[facei];

	Uif[facei] = 0.25*sqrt(8*e*Ti_/pi/mi_)*n[facei] + Ui[faceCelli];
    }

    operator == (Uif);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::ionVelocity::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Ti")
	    << Ti_ << token::END_STATEMENT << nl;
    os.writeKeyword("mi")
	    << mi_ << token::END_STATEMENT << nl;
    os.writeKeyword("variable")
	    << variable_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        ionVelocity
    );
}

// ************************************************************************* //
