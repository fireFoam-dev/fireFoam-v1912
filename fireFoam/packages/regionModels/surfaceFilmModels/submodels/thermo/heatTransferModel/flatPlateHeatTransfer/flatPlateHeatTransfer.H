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
    Foam::flatPlateHeatTransfer

Description
    Heat transfer model based on Al Khalil's model

SourceFiles
    flatPlateHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef flatPlateHeatTransfer_H
#define flatPlateHeatTransfer_H

#include "heatTransferModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class flatPlateHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class flatPlateHeatTransfer
:
    public heatTransferModel
{
private:

    // Private data

    // Private member functions

        //- Disallow default bitwise copy construct
        flatPlateHeatTransfer(const flatPlateHeatTransfer&);

        //- Disallow default bitwise assignment
        void operator=(const flatPlateHeatTransfer&);

protected:
    // Protected data
        //- Length scale / [m]
        scalar L_;

        //- Heat transfer coefficient - film-wall region [W/m2/K]
        volScalarField htcConvFilm_;


    // Protected member functions

        //- Return Nusselt number as a function of Reynolds and Prandtl numbers
        scalar Nu(const scalar Re, const scalar Sc) const;

public:

    //- Runtime type information
    TypeName("flatPlate");


    // Constructors

        //- Construct from surface film model and dictionary
        flatPlateHeatTransfer
        (
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~flatPlateHeatTransfer();


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
