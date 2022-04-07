/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

#include "distributionContactAngleForceFF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(distributionContactAngleForceFF, 0);
addToRunTimeSelectionTable(force, distributionContactAngleForceFF, dictionary);



void distributionContactAngleForceFF::updateTheta()
{
    scalar time=filmModel_.time().value();

    forAll(contactAngle_,celli){
        if(time > timeOld_[celli]+timeInterval_[celli]){
            contactAngleOld_[celli]=contactAngleNew_[celli];
            timeOld_[celli]=timeOld_[celli]+timeInterval_[celli];
            contactAngleNew_[celli]=distribution_->sample();
            timeInterval_[celli]=timeIntervalDistribution_->sample();
        }
        scalar f = (time - timeOld_[celli])/(timeInterval_[celli]);
        contactAngle_[celli] = (1.0-f)*contactAngleOld_[celli]+f*contactAngleNew_[celli];
    }

    contactAngleNew_.correctBoundaryConditions();
    contactAngleOld_.correctBoundaryConditions();
    contactAngle_.correctBoundaryConditions();
    timeInterval_.correctBoundaryConditions();

    return;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distributionContactAngleForceFF::distributionContactAngleForceFF
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    rndGen_(),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    ),
    timeIntervalDistribution_ // kvm
    (
        distributionModel::New
        (
            coeffDict_.subDict("timeIntervalDistribution"),
            rndGen_
        )
    ),
    contactAngle_ // kvm
    (
        IOobject
        (
            typeName + ":contactAngle",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    contactAngleOld_ // kvm
    (
        IOobject
        (
            typeName + ":contactAngleOld",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    contactAngleNew_ // kvm
    (
        IOobject
        (
            typeName + ":contactAngleNew",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    timeOld_ // kvm
    (
        IOobject
        (
            typeName + ":timeOld",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    timeInterval_ // kvm
    (
        IOobject
        (
            typeName + ":timeInterval",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
     )
{
    forAll(contactAngleOld_,celli){ // kvm
        contactAngleOld_[celli] = distribution_->sample();
        timeOld_[celli] = film.time().value();
    }
    forAll(contactAngleNew_,celli){ // kvm
        contactAngleNew_[celli] = distribution_->sample();
        timeInterval_[celli] = timeIntervalDistribution_->sample();
    }
    // const kinematicSingleLayer& film = // kvm
    //     dynamic_cast<const kinematicSingleLayer&>(filmModel_);
    Info << "dx " << sqrt(film.magSf()[0]) << nl;




}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

distributionContactAngleForceFF::~distributionContactAngleForceFF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> distributionContactAngleForceFF::theta()
{

    updateTheta();

    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":theta",
                filmModel_.time().timeName(),
                filmModel_.regionMesh()
            ),
            filmModel_.regionMesh(),
            dimensionedScalar("0", dimless, 0)
        )
    );

    volScalarField& theta = ttheta.ref();
    volScalarField::Internal& thetai = theta.ref();

    forAll(thetai, celli)
    {
        thetai[celli] = contactAngle_[celli];
    }

    forAll(theta.boundaryField(), patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            fvPatchField<scalar>& thetaf = theta.boundaryFieldRef()[patchi];
            const fvPatchField<scalar>& contactAnglef = contactAngle_.boundaryField()[patchi];

            forAll(thetaf, facei)
            {
                thetaf[facei] = contactAnglef[facei];
            }
        }
    }

    return ttheta;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
