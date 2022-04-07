/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::alKhalilHeatTransfer

Description
    Heat transfer model based on Al Khalil's model

SourceFiles
    alKhalilHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef alKhalilHeatTransfer_H
#define alKhalilHeatTransfer_H

#include "heatTransferModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class alKhalilHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class alKhalilHeatTransfer
:
    public heatTransferModel
{
private:

    // Private data

        //- Constant heat transfer coefficient [W/m2/K]
        scalar c0_;

        //- Heat transfer coefficient - film-wall region [W/m2/K]
        volScalarField htcConvFilm_;

    // Private member functions

        //- Disallow default bitwise copy construct
        alKhalilHeatTransfer(const alKhalilHeatTransfer&);

        //- Disallow default bitwise assignment
        void operator=(const alKhalilHeatTransfer&);


public:

    //- Runtime type information
    TypeName("alKhalil");


    // Constructors

        //- Construct from surface film model and dictionary
        alKhalilHeatTransfer
        (
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~alKhalilHeatTransfer();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();

            //- Return the heat transfer coefficient [W/m2/K]
            virtual tmp<volScalarField> h() const;
            //virtual void h() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
