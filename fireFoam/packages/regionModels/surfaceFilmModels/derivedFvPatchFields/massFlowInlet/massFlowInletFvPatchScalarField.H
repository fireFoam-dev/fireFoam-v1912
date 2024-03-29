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
    Foam::massFlowInletFvPatchScalarField

Description
     inlet boundary condition for patches where the film height is
    specified. The inflow  is obtained from the flux with a direction
    normal to the patch faces.

SourceFiles
    massFlowInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef massFlowInletFvPatchScalarField_H
#define massFlowInletFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class massFlowInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class massFlowInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of film region
        word filmRegionName_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of density field
        word muName_;

        //- Name of film height field
        word deltafName_;

        //- gamma
        scalar gamma_;

        //- mean film thickness
        scalar mean_;

        //- amplitude of oscillation
        scalar amplitude_;

        //- spatial frequency
        scalar spatialFrequency_;

        //- time frequency
        scalar temporalFrequency1_;

        //- time frequency
        scalar temporalFrequency2_;


public:

    //- Runtime type information
    TypeName("massFlowInlet");


    // Constructors

        //- Construct from patch and internal field
        massFlowInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        massFlowInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given massFlowInletFvPatchScalarField
        //  onto a new patch
        massFlowInletFvPatchScalarField
        (
            const massFlowInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        massFlowInletFvPatchScalarField
        (
            const massFlowInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new massFlowInletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        massFlowInletFvPatchScalarField
        (
            const massFlowInletFvPatchScalarField&,
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
                new massFlowInletFvPatchScalarField(*this, iF)
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

            //- Return the name of mu
            const word& muName() const
            {
                return muName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& muName()
            {
                return muName_;
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

        virtual void operator=(const fvPatchField<scalar>& pvf);

protected:

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
