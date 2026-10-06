/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::electronTemperature

Description
    Wall condition for the electron temperature, after Lymberopoulos and
    Economou (J. Phys. D 28 (1995) 727, table 1): the electron energy flux
    to the wall is the enthalpy carried by the electrons,

        q.n = 5/2 k Te Gamma_out - 5/2 k (Gamma_se Tse + Gamma_FE TFN),

    with Gamma_out the flux of plasma electrons to the wall, Gamma_se the
    flux of secondary electrons (seec times the ion flux) and Gamma_FE the
    field-emitted flux. The energy equation carries
    5/2 k Te Gamma_net - kappa dTe/dn through the wall face, with
    Gamma_net = Gamma_out - Gamma_se - Gamma_FE the net electron flux of
    driftDiffusionElectronDensity, so the condition on the temperature is

        kappa dTe/dn = -h (Te - Tem),   h = 5/2 k (Gamma_se + Gamma_FE),

    with Tem the mean temperature of the emitted electrons. It is a mixed
    condition with

        valueFraction = h/(h + kappa/delta),  refValue = Tem,  refGrad = 0.

    h is never negative, so the weight lies between 0 and 1 for any size of
    the cell next to the wall; without emission the condition is a zero
    gradient.

    \verbatim
    electrode
    {
        type            electronTemperature;
        seec            0.05;
        Tse             1000;
        Edepend         true;       // accepted, not used
        TFN             1000;
        field_emission  false;
        field_enhancement_factor 55;
        work_function   5;
    }
    \endverbatim

SourceFiles
    electronTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef electronTemperature_H
#define electronTemperature_H

#include "fvPatchFields.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class electronTemperature Declaration
\*---------------------------------------------------------------------------*/


class electronTemperature
:
    public mixedFvPatchField<scalar>
{
    // Private data
    scalar seec_;
    scalar Tse_;

    bool Edepend_;

    scalar TFN_ ;

    bool FE_;

    scalar beta_;

    scalar wf_;


    public:

    //- Runtime type information
    TypeName("electronTemperature");


    // Constructors

        //- Construct from patch and internal field
        electronTemperature
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronTemperature
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inletOutletFvPatchField onto a new patch
        electronTemperature
        (
            const electronTemperature&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronTemperature
        (
            const electronTemperature&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new electronTemperature(*this)
            );
        }

        //- Construct as copy setting internal field reference
        electronTemperature
        (
            const electronTemperature&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new electronTemperature(*this, iF)
            );
        }


    // Member functions


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<scalar>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
