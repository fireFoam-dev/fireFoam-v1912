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

#include "thermoSingleLayerPw.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedWallPolyPatch.H"
#include "mapDistribute.H"
#include "constants.H" //to get pi
#include <stdio.h>

// Sub-models
#include "injectionModel.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

#include "partialWetting.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayerPw, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerPw, mesh);



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayerPw::read()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerPw::read()" << endl;
    }
    if (thermoSingleLayer::read())
    {
        partialWetting_ = coeffs_.lookupOrDefault<Switch>("partialWetting",false);
        if(partialWetting_){
            const dictionary& subdict=coeffs_.subDict("partialWettingCoeffs");
            subdict.lookup("contactAngleCoefficient") >> contactAngleCoefficient_; 
            subdict.lookup("criticalFilmThickness") >> criticalFilmThickness_; 
            subdict.lookup("contactAngleMax") >> contactAngleMax_; 
            subdict.lookup("contactAngleMin") >> contactAngleMin_; 
            subdict.lookup("contactAngleMean") >> contactAngleMean_; 
            subdict.lookup("contactAngleStdDev") >> contactAngleStdDev_; 
            subdict.lookup("hydrophilic") >> hydrophilic_; 
            /*DEBUG(hydrophilic_)*/
            subdict.lookup("dryToggle") >> dryToggle_; 
            subdict.lookup("wetToggle") >> wetToggle_; 
        }

        absorption_ = coeffs_.lookupOrDefault<Switch>("absorption",false);
        
        if (debug)
        {
            Info<<fmtab.c_str()<< "leaving thermoSingleLayerPw::read()" << endl;
            tabSubtract();
        }
        return true;
    }
    else
    {
        return false;
    }
}

void thermoSingleLayerPw::updateSubmodels()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerPw::updateSubmodels()" << endl;
    }
    /*vaporization and separation are added to rhoSp_*/
    thermoSingleLayer::updateSubmodels();
    if(absorption_){
        updateMassAbsorption();
        rhoSp_ +=massAbsorption_/magSf()/time_.deltaT();
    }

    /*partially wetted treatment*/
    updateOmega();
    updateAlpha();
    updateContactLine();

    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerPw::updateSubmodels()" << endl;
        tabSubtract();
    }
}

tmp<fvVectorMatrix> thermoSingleLayerPw::solveMomentum
(
    const volScalarField& pu,
    const volScalarField& pp
)
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerPw::solveMomentum()" << endl;
    }

    updateSurfaceVelocities();

    contactAngleForce_ = +1.0*contactAngleCoefficient_*sigma_*(1-cos(contactAngle_*constant::mathematical::pi/180))/sqrt(magSf())*contactLine_;

    // Momentum
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(deltaRho_, U_)
      + fvm::div(phi_, U_)
     //TODO: have Andy add this term
     //- delta_*fvm::laplacian(mu_,U_)
     ==
      - USp_
      - rhoSp_*U_
      + forces_.correct(U_)
      + turbulence_->Su(U_)
      + contactAngleForce_
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    if (momentumPredictor_)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
              - fvc::interpolate(delta_)
              * (
                    regionMesh().magSf()
                  * (
                        fvc::snGrad(pu, "snGrad(p)")
                      + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
                      + fvc::snGrad(delta_)*fvc::interpolate(pp)
                    )
                  - (fvc::interpolate(rho_*gTan()) & regionMesh().Sf())
                )
            )
        );

        // Remove any patch-normal components of velocity
        U_ -= nHat()*(nHat() & U_);
        U_.correctBoundaryConditions();
    }

    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerPw::solveMomentum()" << endl;
        tabSubtract();
    }
    return tUEqn;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerPw::thermoSingleLayerPw
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    thermoSingleLayer(modelType, mesh, g, regionType),
    partialWetting_(false),
    absorption_(false),
    contactAngleMean_(0.0), 
    contactAngleStdDev_(0.0), 
    contactAngleMax_(180.0), 
    contactAngleMin_(0.0), 
    contactAngleCoefficient_(0.0), 
    contactAngleFromFile_(false), //kvm
    criticalFilmThickness_("criticalFilmThickness",dimLength,0.0), 
    //hydrophilic_(false),
    dryToggle_(0.01),
    wetToggle_(0.5),

    omega_ //omega=0 -> dry; omega=1 -> wet;
    (
         IOobject
         (
             "omega",
             time().timeName(),
             regionMesh(),
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         ),
         regionMesh()
     ),
     contactAngle_
     (
         IOobject
         (
             "contactAngle",
             time().timeName(),
             regionMesh(),
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE
         ),
         regionMesh(),
         0.0,
         zeroGradientFvPatchVectorField::typeName
     ),
     contactAngleForce_
     (
         IOobject
         (
             "contactAngleForce",
             time().timeName(),
             regionMesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         regionMesh(),
         dimensionedVector("zero", dimMass/sqr(dimTime)/dimLength, vector::zero),
         zeroGradientFvPatchVectorField::typeName
     ),
    contactLine_
    (
        IOobject
        (
            "contactLine",
            time_.timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    timeWetted_
    (
        IOobject
        (
            "timeWetted",
            time_.timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTime, 0.0)
    ),
    massAbsorption_
    (
        IOobject
        (
            "massAbsorption",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    cummulativeMassAbsorption_
    (
        IOobject
        (
            "cummulativeMassAbsorption",
            time_.timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerPw::~thermoSingleLayerPw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void thermoSingleLayerPw::info()
{
    thermoSingleLayer::info();
    if(absorption_){
        Info<< indent << "absorbed mass         = "
            << gSum((cummulativeMassAbsorption_*magSf())())/gSum(magSf())<< nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
