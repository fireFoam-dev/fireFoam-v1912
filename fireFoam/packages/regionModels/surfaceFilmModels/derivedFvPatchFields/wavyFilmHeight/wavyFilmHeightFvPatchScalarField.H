/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::wavyFilmHeightFvPatchScalarField

Description
    Film height boundary condition that imposes a sinusoidal profile on top of
    a mean film height

SourceFiles
    wavyFilmHeightFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wavyFilmHeightFvPatchScalarField_H
#define wavyFilmHeightFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class wavyFilmHeightFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wavyFilmHeightFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Mean film height
        autoPtr<TimeFunction1<scalar> > mean_;

        //- Perturbation amplitude [m]
        autoPtr<TimeFunction1<scalar> > a_;

        //- Perturbation frequency [rad/s/m]
        autoPtr<TimeFunction1<scalar> > omega_;


public:

    //- Runtime type information
    TypeName("wavyFilmHeight");


    // Constructors

        //- Construct from patch and internal field
        wavyFilmHeightFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wavyFilmHeightFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wavyFilmHeightFvPatchScalarField
        //  onto a new patch
        wavyFilmHeightFvPatchScalarField
        (
            const wavyFilmHeightFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wavyFilmHeightFvPatchScalarField
        (
            const wavyFilmHeightFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wavyFilmHeightFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wavyFilmHeightFvPatchScalarField
        (
            const wavyFilmHeightFvPatchScalarField&,
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
                new wavyFilmHeightFvPatchScalarField(*this, iF)
            );
        }


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
