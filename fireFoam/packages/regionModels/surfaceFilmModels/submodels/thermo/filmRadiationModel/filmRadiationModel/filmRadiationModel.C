/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmRadiationModel, 0);
defineRunTimeSelectionTable(filmRadiationModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmRadiationModel::filmRadiationModel(surfaceFilmRegionModel& film)
:
    filmSubModelBase(film)
{}


filmRadiationModel::filmRadiationModel
(
    const word& modelType,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmSubModelBase(film, dict, typeName, modelType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

filmRadiationModel::~filmRadiationModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// luwi
tmp<volScalarField> filmRadiationModel::Shs(const volScalarField& qemPyro)
{
    return this->Shs();
}
// luwi
const tmp<volScalarField> filmRadiationModel::transmissivity() const
{
    tmp<volScalarField> tOne
    (
        new volScalarField
        (
            IOobject
            (
                "transmissivityOne",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
	        dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tOne;
}

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
