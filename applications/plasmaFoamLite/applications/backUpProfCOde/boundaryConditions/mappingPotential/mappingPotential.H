/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd
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

Class
    Foam::compressible::
        turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer
    to be used for in multiregion cases with two phase Euler system

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        qrNbr      | name of the radiative flux in the nbr region | no | none
        qr         | name of the radiative flux in this region | no | none
        region     | region to which this BC belongs | yes
        otherPhase | name of the vapour phase in the fluid region | yes
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
    \endtable

    Example of the boundary condition specification on the fluid region:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureTwoPhaseRadCoupledMixed;
        Tnbr            T;
        qrNbr           none;
        qr              none;
        kappaMethod     phaseSystem;
        region          fluid;
        otherPhase      gas;
        value           uniform 300;
    }
    \endverbatim

    Example of the boundary condition specification on the solid region:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureTwoPhaseRadCoupledMixed;
        Tnbr            T.liquid;
        qrNbr           none;
        qr              none;
        kappaMethod     solidThermo;
        region          solid;
        otherPhase      gas;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.


SourceFiles
    turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/
#ifndef mappingPotential_H
#define mappingPotential_H
#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class mappingPotential
:
    public mixedFvPatchScalarField
{
public:

     // Public enumerations
        enum KMethodType
        {
            mtSolidThermo,
            mtLookup,
            mtPhaseSystem
        };

     // Data types
        enum regionType
        {
            solid,
            fluid
        };

private:
          const word nbrFieldName_;
          const word nbrPatchName_;

public:
    TypeName("mappingPotential");


    // Constructors

        mappingPotential
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
        mappingPotential
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        mappingPotential
        (

            const mappingPotential&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mappingPotential
                (
                    *this
                )
            );
        }

        mappingPotential
        (

            const mappingPotential&,
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

                new mappingPotential
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions
        const fvPatchScalarField& sigma() const;
        const fvPatchVectorField& EField() const;
        const fvPatchScalarField& epsilon() const;
        const fvPatchScalarField& phi() const;
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
