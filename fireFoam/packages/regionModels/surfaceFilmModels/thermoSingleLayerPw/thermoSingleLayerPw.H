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
    Foam::thermoSingleLayerPw

Description
    Kinematic form of single-cell layer surface film model
    with added treatment for partial wetting

SourceFiles
    thermoSingleLayerPw.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayerPw_H
#define thermoSingleLayerPw_H

#include "thermoSingleLayer.H"
#include "Random.H" //kvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class injectionModel;

/*---------------------------------------------------------------------------*\
                   Class thermoSingleLayerPw Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayerPw
:
    public thermoSingleLayer
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayerPw(const thermoSingleLayerPw&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayerPw&);


protected:

    // Protected data

        // Model parameters

            Switch partialWetting_;  

            Switch absorption_;  

            scalar contactAngleMean_;  

            scalar contactAngleStdDev_;  

            scalar contactAngleMax_;  

            scalar contactAngleMin_;  

            scalar contactAngleCoefficient_;  

            Switch contactAngleFromFile_;  

            dimensionedScalar criticalFilmThickness_;  
           
            //Switch hydrophilic_;  
            
            scalar dryToggle_;  
            
            scalar wetToggle_;  
            
            // Fields

                //- Surface wettness
                volScalarField omega_;

                //- Contact Angle / [deg]
                volScalarField contactAngle_;

                //- Contact Angle Force [N/m2]
                volVectorField contactAngleForce_;

        // Fields
        
            //- Velocity / [m/s]
            volVectorField contactLine_;

            //- time since surface initially became wetted
            volScalarField timeWetted_;
            
            //- Mass absorbed into solid
            volScalarField massAbsorption_;

            //- Mass absorbed into solid
            volScalarField cummulativeMassAbsorption_;

    // Protected member functions

        //- Initialise the film model - called on construction
        void initialise();

        //- Read control parameters from dictionary
        virtual bool read();

        //-Contact Angle Related Fields
        void computeContactAngle();
        
        //-Wet-dry field
        void updateOmega();

        //-Wet-dry field
        virtual void updateAlpha();

        //- Update the film sub-models
        virtual void updateSubmodels();

        //- Update the mass absorption
        void updateMassAbsorption();

        //- Update the contact line
        void slottedCircle();

        //- Update the contact line
        virtual void updateContactLine();

        // Equations

            //- Solve for film velocity
            virtual tmp<fvVectorMatrix> solveMomentum
            (
                const volScalarField& pu,
                const volScalarField& pp
            );

public:

    //- Runtime type information
    TypeName("thermoSingleLayerPw");


    // Constructors

        //- Construct from components
        thermoSingleLayerPw
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );


    //- Destructor
    virtual ~thermoSingleLayerPw();


    // Member Functions

        // Evolution

        // Access functions
        
            //- Return the film wettness (0=dry, 1=wet)
            virtual const volScalarField& omega() const;

            //- Return the pyrolysis virgin fraction
            // virtual const volScalarField& Yv() const;

            
       // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
