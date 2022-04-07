/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::filmHeightInletVelocityWavyFvPatchVectorField

Description
    VelocityWavy inlet boundary condition for patches where the film height is
    specified. The inflow velocity is obtained from the flux with a direction
    normal to the patch faces.

SourceFiles
    filmHeightInletVelocityWavyFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef filmHeightInletVelocityWavyFvPatchVectorField_H
#define filmHeightInletVelocityWavyFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class filmHeightInletVelocityWavyFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class filmHeightInletVelocityWavyFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of film height field
        word deltafName_;


public:

    //- Runtime type information
    TypeName("filmHeightInletVelocityWavy");


    // Constructors

        //- Construct from patch and internal field
        filmHeightInletVelocityWavyFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmHeightInletVelocityWavyFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given filmHeightInletVelocityWavyFvPatchVectorField
        //  onto a new patch
        filmHeightInletVelocityWavyFvPatchVectorField
        (
            const filmHeightInletVelocityWavyFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        filmHeightInletVelocityWavyFvPatchVectorField
        (
            const filmHeightInletVelocityWavyFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new filmHeightInletVelocityWavyFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        filmHeightInletVelocityWavyFvPatchVectorField
        (
            const filmHeightInletVelocityWavyFvPatchVectorField&,
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
                new filmHeightInletVelocityWavyFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of deltaf
            const word& deltafName() const
            {
                return deltafName_;
            }

            //- Return reference to the name of df to allow adjustment
            word& deltafName()
            {
                return deltafName_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
