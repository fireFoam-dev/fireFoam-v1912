/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "alKhalilHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
/*#include "thermoSingleLayerFM.H" //kvm*/
#include "thermoSingleLayer.H" //kvm

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{
defineTypeNameAndDebug(alKhalilHeatTransfer, 0);
addToRunTimeSelectionTable ( heatTransferModel, alKhalilHeatTransfer, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alKhalilHeatTransfer::alKhalilHeatTransfer
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

alKhalilHeatTransfer::~alKhalilHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void alKhalilHeatTransfer::correct()
{
    /*const thermoSingleLayerFM& film = refCast<const thermoSingleLayerFM>(filmModel_); //kvm*/
    const thermoSingleLayer& film = refCast<const thermoSingleLayer>(filmModel_); //kvm
    
    // retrieve fields from film model
    const scalarField& delta = film.delta();
    //const scalarField& T = film.T();
    //const scalarField& Tw = film.Tw();
    const scalarField& rho = film.rho();
    const scalarField& mu = film.mu();
    //const scalarField& magSf = film.magSf();
    const scalarField kappa = film.kappa();
    const vectorField dU = film.Uw() - film.Us();
    
    forAll(htcConvFilm_, cellI)
    {
        // Film region density [kg/m3]
        const scalar rhoc = rho[cellI];
    
        // Film region viscosity [Pa.s]
        const scalar muc = mu[cellI];
    
        // Film region thickness [m]
        const scalar deltac = delta[cellI];
    
        // Reynolds number
        const scalar Rec = rhoc*mag(dU[cellI])*deltac/muc;
    
        // thermal conductivity [W/m/K]
        const scalar kappac = kappa[cellI];
    
        // Nusselt number
        //const scalar Nu = 2.63+0.000143*Rec; //constant Tw
        const scalar Nu = 3.20+0.000237*Rec; //constant q"
    
        // heat transfer coefficient [W/m2/K]
        const scalar hm = Nu*kappac/(deltac+ROOTVSMALL);
    
    //    htcConvFilm_[cellI]=min(hm,5e4);
        htcConvFilm_[cellI]=min(hm,1e4);//based on minimum film thickness of 0.0002, k=0.6, NU=3.1, h = Nu k / delta
    }
}


tmp<volScalarField> alKhalilHeatTransfer::h() const
{
    return htcConvFilm_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
