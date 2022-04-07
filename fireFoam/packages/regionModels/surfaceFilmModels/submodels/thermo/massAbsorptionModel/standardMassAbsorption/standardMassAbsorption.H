/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::standardMassAbsorption

Description
    Standard phase change model with modification for boiling

SourceFiles
    standardMassAbsorption.C

\*---------------------------------------------------------------------------*/

#ifndef standardMassAbsorption_H
#define standardMassAbsorption_H

#include "massAbsorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class standardMassAbsorption Declaration
\*---------------------------------------------------------------------------*/

class standardMassAbsorption
:
    public massAbsorptionModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        standardMassAbsorption(const standardMassAbsorption&);

        //- Disallow default bitwise assignment
        void operator=(const standardMassAbsorption&);


protected:

    // Protected data

        //- Minimum film height for model to be active
        const scalar deltaMin_;


    // Protected member functions

        //- Mass absorbed into solid
        volScalarField massAbsorption_;

        //- Mass absorbed into solid
        volScalarField cummulativeMassAbsorption_;



public:

    //- Runtime type information
    TypeName("standardMassAbsorption");


    // Constructors

        //- Construct from surface film model
        standardMassAbsorption
        (
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~standardMassAbsorption();


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
