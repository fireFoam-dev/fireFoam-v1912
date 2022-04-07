/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::bubblingGasification

Description
    Standard phase change model with modification for boiling

SourceFiles
    bubblingGasification.C

\*---------------------------------------------------------------------------*/

#ifndef bubblingGasification_H
#define bubblingGasification_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class bubblingGasification Declaration
\*---------------------------------------------------------------------------*/

class bubblingGasification
:
    public phaseChangeModel
{
protected:

    // Protected data

        //- Minimum film height for model to be active
        const scalar deltaMin_;

        //- Length scale / [m]
        const scalar L_;

        //- Boiling temperature factor / []
        //  Used to set max limit on temperature to Tb*TbFactor
        const scalar TbFactor_;

        const scalar scaling_; // kvm

        //- Arrhenius rates...
        const scalar A_;  // AK pre-exponential (1/s) for single step global reaction
        const scalar Ta_; // AK activation temperature (K)
        const scalar b_;  // AK temperature dependence exponent


        scalar dMassMax_; // kvm
        // scalarField dMassMax_;

    // Protected member functions

    //- No copy construct
        bubblingGasification(const bubblingGasification&) = delete;

        //- No copy assignment
        void operator=(const bubblingGasification&) = delete;

public:

    //- Runtime type information
    TypeName("bubblingGasification");


    // Constructors

        //- Construct from surface film model
        bubblingGasification
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~bubblingGasification() = default;


    // Member Functions

        // Evolution

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
