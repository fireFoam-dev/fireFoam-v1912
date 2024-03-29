/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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
    Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField

Description
    Describes a volumetric/mass flow normal vector boundary condition by its
    magnitude as an integral over its area.

    The inlet mass flux is taken from the neighbor region.

    phi is used to determine if the flow is compressible or incompressible.

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the velocity when applying the
    mass basis.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            flowRateInletVelocityPyrolysisCoupled;
        phi             phi;
        rho             rho;
        neigPhi         neigPhiName_;    // Volumetric/mass flow rate [m3/s or kg/s]
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    flowRateInletVelocityPyrolysisCoupledFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletVelocityPyrolysisCoupledFvPatchVectorField_H
#define flowRateInletVelocityPyrolysisCoupledFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class flowRateInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the neighbor flux setting the inlet mass flux
        word nbrPhiName_;

        //- Name of the local mass flux
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Heat of combustion of pyrolysate
        const scalar hocSolid_;



public:

   //- Runtime type information
   TypeName("flowRateInletVelocityPyrolysisCoupled");


   // Constructors

        //- Construct from patch and internal field
        flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
        //  onto a new patch
        flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
        (
            const flowRateInletVelocityPyrolysisCoupledFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
        (
            const flowRateInletVelocityPyrolysisCoupledFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletVelocityPyrolysisCoupledFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
        (
            const flowRateInletVelocityPyrolysisCoupledFvPatchVectorField&,
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
                new flowRateInletVelocityPyrolysisCoupledFvPatchVectorField(*this, iF)
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
