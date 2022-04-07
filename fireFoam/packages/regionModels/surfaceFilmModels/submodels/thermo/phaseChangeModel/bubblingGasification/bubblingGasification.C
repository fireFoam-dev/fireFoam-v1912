/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "bubblingGasification.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoSingleLayer.H"
#include "zeroField.H"
#include "constants.H" // ak

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(bubblingGasification, 0);

addToRunTimeSelectionTable
(
    phaseChangeModel,
    bubblingGasification,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bubblingGasification::bubblingGasification
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    phaseChangeModel(typeName, film, dict),
    deltaMin_(coeffDict_.get<scalar>("deltaMin")),
    L_(coeffDict_.get<scalar>("L")),
    TbFactor_(coeffDict_.lookupOrDefault<scalar>("TbFactor", 1.1)),
    scaling_(coeffDict_.lookupOrDefault<scalar>("scaling",1.0)), // kvm
    A_(coeffDict_.get<scalar>("A")), // AK
    Ta_(coeffDict_.get<scalar>("Ta")), // AK
    b_(coeffDict_.lookupOrDefault<scalar>("b",1.0)), // AK
    dMassMax_(coeffDict_.lookupOrDefault<scalar>("dMassMax",5e-6)) // kvm
{ // kvm
    Info << "Mass transfer convective scaling set to " << scaling_ << nl; // kvm
} // kvm



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void bubblingGasification::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    // Set local thermo properties
    const SLGThermo& thermo = film.thermo();
    const filmThermoModel& filmThermo = film.filmThermo();
    const label vapId = thermo.carrierId(filmThermo.name());

    // Retrieve fields from film model
    const scalarField& delta = film.delta();
    const scalarField& pInf = film.pPrimary();
    const scalarField& T = film.T();
    const scalarField& kappa = film.kappa(); // kvm
    const scalarField& rho = film.rho();
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& magSf = film.magSf();
    // Info << "uprimary " << film.UPrimary();
    // Info << "us " << film.Us();
    const scalarField limMass
    (
        max(scalar(0), availableMass - deltaMin_*rho*magSf)
    );
    // const scalarField qRad(film.qRad());
    
    static scalarField dMassPrev(dMass.size(),0.0); // kvm

    scalar boilingMass = 0.0; // kvm
    scalar evaporationMass = 0.0; // kvm

    scalar dMassMax = -GREAT; // kvm
    scalar dMassMin =  GREAT; // kvm

    forAll(dMass, celli)
    {
        if (delta[celli] > deltaMin_)
        {
            // Cell pressure [Pa]
            const scalar pc = max(min(pInf[celli],101325*1.1),101325*0.9); // kvm

            // Calculate the boiling temperature
            const scalar Tb = filmThermo.Tb(pc);

            // estimate surface temperature based on energy balance
            const scalar Tf = T[celli]; // kvm

            // Local temperature - impose lower limit of 200 K for stability
            const scalar Tloc = min(TbFactor_*Tb, max(200.0, Tf)); // kvm

            // Latent heat [J/kg]
            const scalar hVap = filmThermo.hl(pc, Tloc);

            // Calculate mass transfer

            // Use Arrhenius rates here
            const scalar A(A_); // 1/s
            const scalar Ta(Ta_); // K
 
            const scalar omega = A*exp(-Ta/Tf); // 1/s
 
            dMass[celli] = delta[celli]*omega*dt*rho[celli]*magSf[celli];

            dMass[celli] = min(limMass[celli], max(0.0, dMass[celli]));
            // try to under-relax the vaporization rate
            dMass[celli] = min(dMassPrev[celli]*1.05+SMALL, max(0.0, dMass[celli])); // kvm
            // give a hard limit on vaporization rate, 
            // TODO: ideally should be based on cell volume and time step
            // dMass[celli] = min(5e-5, dMass[celli]);
            //dMass[celli] = min(dMassMax_, dMass[celli]); // kvm

            boilingMass += dMass[celli];

            // dMass[celli] = min(5e-7, dMass[celli]);
            dMassMax = max(dMassMax,dMass[celli]); // kvm
            dMassMin = min(dMassMin,dMass[celli]); // kvm

            // Info << "dMassa " << dMass[celli] << endl;
    
            dEnergy[celli] = dMass[celli]*hVap;
            dMassPrev[celli] = dMass[celli]; // kvm
        }
    }
    reduce(dMassMax,maxOp<scalar>()); // kvm
    // dMassMax_ = max(dMassMax,dMassMax_);
    reduce(dMassMin,minOp<scalar>()); // kvm
    Info<<"bubblingGasification::dMassMax" << tab << film.time().timeName() << tab << dMassMax <<endl; // kvm
    // Info<<"bubblingGasification::dMassMax_" << tab << film.time().timeName() << tab << dMassMax_ <<endl; // kvm
    Info<<"bubblingGasification::dMassMin" << tab << film.time().timeName() << tab << dMassMin <<endl; // kvm
    reduce(boilingMass,sumOp<scalar>()); // kvm
    Info << "bubblingGasification::massLossTotal " << boilingMass << nl; // kvm
}


//void bubblingGasification::correctModel
//(
//    const scalar dt,
//    scalarField& availableMass,
//    scalarField& dMass,
//    scalarField& dEnergy
//)
//{
//    if (YInfZero_)
//    {
//        correctModel(dt, availableMass, dMass, dEnergy, zeroField());
//    }
//    else
//    {
//        const thermoSingleLayer& film = filmType<thermoSingleLayer>();
//        const label vapId = film.thermo().carrierId(film.filmThermo().name());
//        const scalarField& YInf = film.YPrimary()[vapId];
//
//        correctModel(dt, availableMass, dMass, dEnergy, YInf);
//    }
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
