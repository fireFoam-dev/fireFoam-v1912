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

#include "massAbsorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(massAbsorptionModel, 0);
defineRunTimeSelectionTable(massAbsorptionModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

massAbsorptionModel::massAbsorptionModel
(
    surfaceFilmRegionModel& film
)
:
    filmSubModelBase(film),
    latestMassAbs_(0.0),
    totalMassAbs_(0.0)
{}


massAbsorptionModel::massAbsorptionModel
(
    const word& modelType,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmSubModelBase(film, dict, typeName, modelType),
    latestMassAbs_(0.0),
    totalMassAbs_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

massAbsorptionModel::~massAbsorptionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void massAbsorptionModel::correct
(
    const scalar dt,
    scalarField& availableMass,
    volScalarField& dMass,
    volScalarField& dEnergy
)
{
    if (!active())
    {
        return;
    }

    correctModel
    (
        dt,
        availableMass,
        dMass,
        dEnergy
    );

    latestMassAbs_ = sum(dMass.primitiveField());
    totalMassAbs_ += latestMassAbs_;

    availableMass -= dMass;
    dMass.correctBoundaryConditions();
}


void massAbsorptionModel::info(Ostream& os) const
{
    const scalar massAbsRate =
        returnReduce(latestMassAbs_, sumOp<scalar>())
       /filmModel_.time().deltaTValue();

    os  << indent << "mass absorption    = " 
        << returnReduce(totalMassAbs_, sumOp<scalar>()) << nl
        << indent << "absorption rate    = " << massAbsRate << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam

// ************************************************************************* //
