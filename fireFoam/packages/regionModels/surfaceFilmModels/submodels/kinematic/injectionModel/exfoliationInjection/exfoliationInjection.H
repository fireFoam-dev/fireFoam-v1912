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
    Foam::exfoliationInjection

Description
    Film Exfoliation mass transfer model.

    If the film encounters a non-zero value of exfoliation, 
    then a certain percentage of the mass is expelled as 
    droplets. 
    If the mass exceeds that needed to generate a valid parcel, the
    equivalent mass is removed from the film.

    New parcel diameters are sampled from a PDF.

SourceFiles
    exfoliationInjection.C

\*---------------------------------------------------------------------------*/

#ifndef exfoliationInjection_H
#define exfoliationInjection_H

#include "injectionModel.H"
#include "distributionModel.H"
#include "cachedRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class exfoliationInjection Declaration
\*---------------------------------------------------------------------------*/

class exfoliationInjection
:
    public injectionModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        exfoliationInjection(const exfoliationInjection&);

        //- Disallow default bitwise assignment
        void operator=(const exfoliationInjection&);


protected:

    // Protected data

        //- Stable film thickness - drips only formed if thickness
        //  execeeds this threhold value
        scalar deltaStable_;

        // - Fraction of available mass to be expelled as droplets
        scalar expelFraction_;

        // - At what value of Yv to expel surface film
        scalar expelThreshold_;

        //- Number of particles per parcel
        scalar particlesPerParcel_;

        //- Random number generator
        cachedRandom rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModels::distributionModel>
            parcelDistribution_;

        //- Diameters of particles to inject into the exfoliation
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("exfoliationInjection");


    // Constructors

        //- Construct from surface film model
        exfoliationInjection
        (
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~exfoliationInjection();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToInject,
                scalarField& diameterToInject
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
