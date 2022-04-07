/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    thermoSingleLayerFmBox

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15 K

SourceFiles
    thermoSingleLayerFmBox.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayerFmBox_H
#define thermoSingleLayerFmBox_H

#include "thermoSingleLayerPw.H"
#include "SLGThermo.H"

#define DEBUG(x) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
#define TRACE(s) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class heatTransferModel;
class phaseChangeModel;

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayerFmBox Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayerFmBox
:
    public thermoSingleLayerPw
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayerFmBox(const thermoSingleLayerFmBox&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayerFmBox&);


protected:

    // Protected data

        //- Sum of trapped mass per tier
        DynamicList< scalar > trapSum;

        //- Activate trapped mass model
        Switch trappedMassModel_;

        //- Number of boxes in array
        label  nBoxes_;

        //- Number of tiers in array
        label  nTiers_;

        //- Empirical parameter
        scalar trapFraction_;

        //- Empirical parameter
        scalar trapProfileFactor_;

        //- Maximum mass trapped per box
        dimensionedScalar trapMax_;

        //- Mass lost due to entrapment
        volScalarField trappedMass_;

        //- Total mass lost
        volScalarField massLost_;

        //- Mass absorbed into solid
        volScalarField rhoSpLost_;

        //- Box
        volScalarField box_;

        //- Tier
        volScalarField tier_;

        //- Face
        volScalarField side_;

        //- Top Edge
        volScalarField topEdge_;
        //
        //- Bottom Edge
        volScalarField bottomEdge_;

        // file for film extent output
        //OFstream outFilmExtent_;

    // Protected member functions

        //- Initialise the film model - called on construction
        void initialise();

        //- Read control parameters from dictionary
        virtual bool read();

        //-Contact Angle Related Fields
        void assignBoxTierSide();//kvm
        
        //-Contact Angle Related Fields
        void integrateMass();//kvm

        //- Update the film sub-models
        virtual void updateSubmodels();

        //- Update the mass absorption
        void updateMassLost();

public:

    //- Runtime type information
        TypeName("thermoSingleLayerFmBox");

    // Constructors

        //- Construct from components
        thermoSingleLayerFmBox
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );


    //- Destructor
        virtual ~thermoSingleLayerFmBox();

    // I-O

        //- Provide some feedback
        virtual void info();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
