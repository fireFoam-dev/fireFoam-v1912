/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::filmDiversionForce

Description
    Film diversion force

SourceFiles
    filmDiversionForce.C

\*---------------------------------------------------------------------------*/

#ifndef filmDiversionForce_H
#define filmDiversionForce_H

#include "force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                  Class filmDiversionForce Declaration
\*---------------------------------------------------------------------------*/

class filmDiversionForce
:
    public force
{
private:

        //- film diversion coefficient
        scalar coeff_;

    // Private member functions

        //- No copy construct
        filmDiversionForce(const filmDiversionForce&) = delete;

        //- No copy assignment
        void operator=(const filmDiversionForce&) = delete;


public:

    //- Runtime type information
    TypeName("filmDiversion");


    // Constructors

        //- Construct from surface film model
        filmDiversionForce
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmDiversionForce();


    // Member Functions

        // Evolution

            //- Correct
            virtual tmp<fvVectorMatrix> correct(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
