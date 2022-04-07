/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "standardMassAbsorption.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoSingleLayer.H"
#include "specie.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardMassAbsorption, 0);

addToRunTimeSelectionTable
(
    massAbsorptionModel,
    standardMassAbsorption,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardMassAbsorption::standardMassAbsorption
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    massAbsorptionModel(typeName, film, dict),
    deltaMin_(readScalar(coeffDict_.lookup("deltaMin"))),
    massAbsorption_
    (
        IOobject
        (
            "massAbsorption",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    cummulativeMassAbsorption_
    (
        IOobject
        (
            "cummulativeMassAbsorption",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    )

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardMassAbsorption::~standardMassAbsorption()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardMassAbsorption::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    if (debug)
    {
        Info<< "standardMassAbsorption::correctModel()" << endl;
    }

    const thermoSingleLayer& film = refCast<const thermoSingleLayer>(filmModel_);

    // retrieve fields from film model
    const scalarField& rho = film.rho();
    const scalarField& magSf = film.magSf();
    const scalarField limMass
    (
        max(scalar(0.0), availableMass - deltaMin_*rho*magSf)
    );
    const scalarField& alpha = film.alpha();
    const scalarField& T = film.T();
    const scalarField& hs = film.hs();

    const scalar T1=20.0+273.15; //K
    const scalar A1=0.0207; //kg/m2
    const scalar n1=0.456;
    const scalar T2=43.0+273.15; //K
    const scalar A2=0.0317; //kg/m2
    const scalar n2=0.465;

    forAll(massAbsorption_,cellI){
        if(alpha[cellI]==1){
            /*linear interpolation for absorption pre-factor*/
            scalar Tloc=max(T[cellI],T1);
            Tloc=min(Tloc,T2);
            scalar A = (Tloc-T1)/(T2-T1)*(A2-A1)+A1;
            scalar n = (Tloc-T1)/(T2-T1)*(n2-n1)+n1;
            const scalar t=max(pow(cummulativeMassAbsorption_[cellI]/A,1/n)*60.0,dt);
            const scalar to=max(SMALL,t-dt);
            massAbsorption_[cellI]=A*(pow(t/60.0,n)-pow((to)/60.0,n))*magSf[cellI];
            cummulativeMassAbsorption_[cellI]+=massAbsorption_[cellI]/magSf[cellI];

            dMass[cellI] = min(limMass[cellI], max(0.0, massAbsorption_[cellI]));
            dEnergy[cellI] = dMass[cellI]*hs[cellI];
        }
        else{
            dMass[cellI]=0.0;
            dEnergy[cellI] = 0.0;
        }
    }
    if (debug)
    {
        Info<< "leaving standardMassAbsorption::correctModel()" << endl;
        
    }




}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
