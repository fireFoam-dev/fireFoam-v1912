/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::nutPyrolysisWallFunctionFvPatchScalarField

Description
    Buoyant Thermal wall function for turbulent kinematic viscosity condition. 
    Based on "Asymptotic analysis of the near-wall region of turbulence natural
    convection flows"
    Holling and Herwig. J. Fluid Mech (2005), vol 541, pp 383-397

SourceFiles
    nutPyrolysisWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutPyrolysisWallFunctionFvPatchScalarField_H
#define nutPyrolysisWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nutPyrolysisWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutPyrolysisWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

        bool correctInternal_;
        scalar Cb_;
        scalar fuelConversionRatio_;

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate yPlus
        //virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;


public:

    //- Runtime type information
    TypeName("nutPyrolysisWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutPyrolysisWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutPyrolysisWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutPyrolysisWallFunctionFvPatchScalarField
        //  onto a new patch
        nutPyrolysisWallFunctionFvPatchScalarField
        (
            const nutPyrolysisWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutPyrolysisWallFunctionFvPatchScalarField
        (
            const nutPyrolysisWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutPyrolysisWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutPyrolysisWallFunctionFvPatchScalarField
        (
            const nutPyrolysisWallFunctionFvPatchScalarField&,
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
                new nutPyrolysisWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;
        
        //- Calculate Uplus
        //scalar UYplus(scalar,scalar,scalar);
        //tmp<scalarField> UYplus() const;
        //tmp<scalarField> UYplus(const scalarField& magUp, const scalarField& y, const scalarField& mf) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
