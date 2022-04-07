/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a wide-band, diffusive radiation
    condition, where the patch temperature is specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        T            | temperature field name  | no          | T
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            wideBandDiffusiveRadiationFF;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField
    Foam::radiationCoupledBaseFF

SourceFiles
    wideBandDiffusiveRadiationFFMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wideBandDiffusiveRadiationFFMixedFvPatchScalarField_H
#define wideBandDiffusiveRadiationFFMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "radiationCoupledBaseFF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
      Class wideBandDiffusiveRadiationFFMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wideBandDiffusiveRadiationFFMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public radiationCoupledBaseFF
{
    // Private data

        //- Name of temperature field
        word TName_;


public:

    //- Runtime type information
    TypeName("wideBandDiffusiveRadiationFF");


    // Constructors

        //- Construct from patch and internal field
        wideBandDiffusiveRadiationFFMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wideBandDiffusiveRadiationFFMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given GreyDiffusiveRadiationMixedFvPatchField
        //  onto a new patch
        wideBandDiffusiveRadiationFFMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationFFMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wideBandDiffusiveRadiationFFMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationFFMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wideBandDiffusiveRadiationFFMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wideBandDiffusiveRadiationFFMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationFFMixedFvPatchScalarField&,
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
                new wideBandDiffusiveRadiationFFMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the temperature field name
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the temperature field name to allow
            //  adjustment
            word& TName()
            {
                return TName_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
