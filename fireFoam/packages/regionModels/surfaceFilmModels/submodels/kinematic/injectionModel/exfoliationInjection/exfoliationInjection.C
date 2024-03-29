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

\*---------------------------------------------------------------------------*/

#include "exfoliationInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(exfoliationInjection, 0);
addToRunTimeSelectionTable(injectionModel, exfoliationInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

exfoliationInjection::exfoliationInjection
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    injectionModel(type(), owner, dict),
    expelFraction_(readScalar(coeffs_.lookup("expelFraction"))),
    expelThreshold_(readScalar(coeffs_.lookup("expelThreshold"))),
    particlesPerParcel_(readScalar(coeffs_.lookup("particlesPerParcel"))),
    rndGen_(),
    parcelDistribution_
    (
        distributionModels::distributionModel::New
        (
            coeffs_.subDict("parcelDistribution"),
            rndGen_
        )
    ),
    diameter_(owner.regionMesh().nCells(), 0.0)
{
    forAll(diameter_, faceI)
    {
        diameter_[faceI] = parcelDistribution_->sample();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

exfoliationInjection::~exfoliationInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void exfoliationInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(this->owner());

    const scalar pi = constant::mathematical::pi;

    // calculate available exfoliation mass
//    tmp<volScalarField> tgNorm(film.gNorm());
//    const scalarField& gNorm = tgNorm();
    const scalarField& magSf = film.magSf();

    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& Yv = film.Yv();

//    const scalarField& arbitraryInjection = film.arbitraryInjection();

    scalarField massDrip(film.regionMesh().nCells(), 0.0);

    forAll(Yv, i)
    {
        if (Yv[i] < expelThreshold_)
        {
            const scalar ddelta = max(0.0, delta[i]);
            massDrip[i] +=
                expelFraction_*min(availableMass[i], max(0.0, ddelta*rho[i]*magSf[i]));
        }
    }


    // Collect the data to be transferred
    forAll(massToInject, cellI)
    {
        scalar rhoc = rho[cellI];
        scalar diam = diameter_[cellI];
        scalar minMass = particlesPerParcel_*rhoc*pi/6*pow3(diam);

        if (massDrip[cellI] > minMass)
        {
            // All drip mass can be injected
            massToInject[cellI] += massDrip[cellI];
            availableMass[cellI] -= massDrip[cellI];

            // Set particle diameter
            diameterToInject[cellI] = diameter_[cellI];

            // Retrieve new particle diameter sample
            diameter_[cellI] = parcelDistribution_->sample();
        }
        else
        {
            // Mass below minimum threshold - cannot be injected
            massToInject[cellI] = 0.0;
            diameterToInject[cellI] = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
