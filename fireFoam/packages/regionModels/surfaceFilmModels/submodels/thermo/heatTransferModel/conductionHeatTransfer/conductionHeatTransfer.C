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

\*---------------------------------------------------------------------------*/

#include "conductionHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "thermoSingleLayer.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{
defineTypeNameAndDebug(conductionHeatTransfer, 0);
addToRunTimeSelectionTable ( heatTransferModel, conductionHeatTransfer, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

conductionHeatTransfer::conductionHeatTransfer
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    heatTransferModel(typeName, film, dict),
    c0_(readScalar(coeffDict_.lookup("c0"))),
    htcConvFilm_
    (
        IOobject
        (
            "htcConvWall",
            filmModel_.time().timeName(),
            filmModel_.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmModel_.regionMesh(),
        dimensionedScalar("zero", dimMass/pow3(dimTime)/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

conductionHeatTransfer::~conductionHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void conductionHeatTransfer::correct()
{
const thermoSingleLayer& film = refCast<const thermoSingleLayer>(filmModel_); //kvm

// retrieve fields from film model
const scalarField& delta = film.delta();
const scalarField kappa = film.kappa();

forAll(htcConvFilm_, cellI)
{
    // Film region thickness [m]
    const scalar deltac = delta[cellI];

    // thermal conductivity [W/m/K]
    const scalar kappac = kappa[cellI];

    // heat transfer coefficient [W/m2/K]
    const scalar hm = 2.0*kappac/(deltac+ROOTVSMALL);

    // based on minimum film thickness of 0.0001, k=0.6, h = 2 * k / delta
    htcConvFilm_[cellI]=min(hm,1.2e4);
}
}


tmp<volScalarField> conductionHeatTransfer::h() const
{
    return htcConvFilm_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
