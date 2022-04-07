/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::distributionContactAngleForceFF

Description
    PDF distribution based film contact angle force

See also
    Foam::regionModels::surfaceFilmModels::contactAngleForceFF
    Foam::distributionModel

SourceFiles
    distributionContactAngleForceFF.C

\*---------------------------------------------------------------------------*/

#ifndef distributionContactAngleForceFF_H
#define distributionContactAngleForceFF_H

#include "contactAngleForce.H"
#include "distributionModel.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                Class distributionContactAngleForceFF Declaration
\*---------------------------------------------------------------------------*/

class distributionContactAngleForceFF
:
    public contactAngleForce
{
    // Private Data

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> distribution_;

        const autoPtr<distributionModel> timeIntervalDistribution_; // kvm

        //- Static contact angle
        volScalarField contactAngle_; // kvm
        volScalarField contactAngleOld_; // kvm
        volScalarField contactAngleNew_; // kvm
        volScalarField timeOld_; // kvm
        volScalarField timeInterval_; // kvm

    // Private member functions

        //- Calculate theta
        void updateTheta();

        //- Disallow default bitwise copy construct
        distributionContactAngleForceFF(const distributionContactAngleForceFF&);

        //- Disallow default bitwise assignment
        void operator=(const distributionContactAngleForceFF&);


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta();


public:

    //- Runtime type information
    TypeName("distributionContactAngleFF");


    // Constructors

        //- Construct from surface film model
        distributionContactAngleForceFF
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distributionContactAngleForceFF();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
