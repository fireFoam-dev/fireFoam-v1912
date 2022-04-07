/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::rampingRadiation

Description
    Film constant radiation model.  The constant radiative flux is specified
    by the user, and operated over a time interval defined by a start time and
    duration.  In addition, a mask can be applied to shield the film from the
    radiation.

SourceFiles
    rampingRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef rampingRadiation_H
#define rampingRadiation_H

#include "filmRadiationModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class rampingRadiation Declaration
\*---------------------------------------------------------------------------*/

class rampingRadiation
:
    public filmRadiationModel
{
private:

    // Private data

        //- Constant radiative flux [kg/s3]
        volScalarField qrConst_;

        //- Radiation mask
        volScalarField mask_;

        //- Absorptivity
        scalar absorptivity_;

        //- 
        // scalar initialValue_;

        //- 
        scalar rampTimeInterval_;

        //- 
        scalar rampStep_;

        //- 
        scalar rampStartTime_;

        //- Incident radiation (used for turbulentTemp.. BC)
        volScalarField qin_;

        //- Time start [s]
        const scalar timeStart_;

        //- Duration [s]
        const scalar duration_;


    // Private member functions

        //- Disallow default bitwise copy construct
        rampingRadiation(const rampingRadiation&);

        //- Disallow default bitwise assignment
        void operator=(const rampingRadiation&);


public:

    //- Runtime type information
    TypeName("rampingRadiation");


    // Constructors

        //- Construct from surface film model and dictionary
        rampingRadiation
        (
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rampingRadiation();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();

            //- Return the radiation sensible enthalpy source
            virtual tmp<volScalarField> Shs();

            //- Return the radiation sensible enthalpy source
            virtual tmp<volScalarField> ShsConst() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
