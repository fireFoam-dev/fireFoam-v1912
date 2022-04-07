/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::noMassAbsorption

Description
    Dummy phase change model for 'none'

SourceFiles
    noMassAbsorption.C

\*---------------------------------------------------------------------------*/

#ifndef noMassAbsorption_H
#define noMassAbsorption_H

#include "massAbsorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class noMassAbsorption Declaration
\*---------------------------------------------------------------------------*/

class noMassAbsorption
:
    public massAbsorptionModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        noMassAbsorption(const noMassAbsorption&);

        //- Disallow default bitwise assignment
        void operator=(const noMassAbsorption&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from surface film model
        noMassAbsorption(surfaceFilmRegionModel& owner, const dictionary& dict);


    //- Destructor
    virtual ~noMassAbsorption();


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
