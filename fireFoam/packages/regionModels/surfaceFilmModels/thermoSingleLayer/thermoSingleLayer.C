/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd
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

#include "thermoSingleLayer.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcFlux.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "massAbsorptionModel.H"
#include "filmRadiationModel.H"
#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, thermoSingleLayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList thermoSingleLayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchi)
    {
        if
        (
            T_.boundaryField()[patchi].fixesValue()
         || bTypes[patchi] == mappedFieldFvPatchField<scalar>::typeName
        )
        {
            bTypes[patchi] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar(hsSp_.dimensions(), Zero);
}


void thermoSingleLayer::correctThermoFields()
{
    // limit temperature to reasonable bounds
    forAll(T_, celli) // kvm
    {
        T_[celli] = max(min(T_[celli],Tmax_),Tmin_); // kvm
    }
    
    rho_ == filmThermo_->rho();
    mu_ == filmThermo_->mu(); // OpenFOAM version missing this, bug, kvm
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();

    forAll(rho_, celli) // kvm
    {
        const scalar Ts = max(min(Ts_[celli],Tmax_),Tmin_); // kvm
        const scalar p = pPrimary_[celli]; // kvm
    }

    rho_.correctBoundaryConditions(); // kvm
    mu_.correctBoundaryConditions(); // kvm
    sigma_.correctBoundaryConditions(); // kvm
    Cp_.correctBoundaryConditions(); // kvm
    kappa_.correctBoundaryConditions(); // kvm

}


void thermoSingleLayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    volScalarField::Boundary& hsBf = hs_.boundaryFieldRef();

    forAll(hsBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchi];
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            hsBf[patchi] == hs(Tp, patchi);
        }
    }
}


void thermoSingleLayer::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // pull char fraction from pyrolysis region
    charFrac_.correctBoundaryConditions();
    //luwi-- pull radiative heat flux from pyrolysis region
    qemPyro_.correctBoundaryConditions();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        if(pyrCoupled_){
        // deprecated    // get pyrolysis internal temperature 
        // deprecated    // and put it on coupled patch boundary of film T_ 

        // deprecated     // UIndirectList<scalar>(Tw_, pp.faceCells()) =
        // deprecated     // 		//T_.boundaryFieldRef()[patchi];
        // deprecated     // 		pyrTemperaturePtr_->boundaryFieldRef()[patchi];
        // deprecated     // TODO: how do I do this with the new coupling?

        // deprecated     scalarList Tpyr(pp.faceCells().size(), 0.0);


        // deprecated     typedef regionModels::pyrolysisModels::pyrolysisModel
        // deprecated         pyrolysisModelType;

        // deprecated     const regionModels::regionModel& pyrolysisRegion =
        // deprecated     db().time().lookupObject<regionModels::regionModel>
        // deprecated         (
        // deprecated             "pyrolysisProperties"
        // deprecated             );

        // deprecated     const pyrolysisModelType& pyrolysisModel =
        // deprecated         dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
        // deprecated     
        // deprecated     pyrolysisModelType& pyrolysis =
        // deprecated         const_cast<pyrolysisModelType&>(pyrolysisModel);
        // deprecated     
        // deprecated     // internal cell tempertaure must be used for stability
        // deprecated     Tpyr = 
        // deprecated         mapRegionPatchInternalField<scalar>
        // deprecated         (
        // deprecated             pyrolysis,
        // deprecated             "T",
        // deprecated             patchi,
        // deprecated             true
        // deprecated             );
        // deprecated     
        // deprecated     UIndirectList<scalar>(Tw_,pp.faceCells()) = 
        // deprecated         Tpyr;

        // deprecated     // get char fraction from pyrolysis model
        // deprecated     scalarList charFrac(pp.faceCells().size(), 0.0);
        // deprecated     charFrac = 
        // deprecated         mapRegionPatchInternalField<scalar>
        // deprecated         (
        // deprecated             pyrolysis,
        // deprecated             "charFrac",
        // deprecated             patchi,
        // deprecated             true
        // deprecated             );

        // deprecated     UIndirectList<scalar>(charFrac_,pp.faceCells()) = 
        // deprecated         charFrac;

        // deprecated     
        // deprecated }
        // else if(1){
        //     // compute Tw based on 0D lumped capacitance model
        //     forAll(Tw_,i){
        //         Tw_[i]=qFilmToWall_[i]*time_.deltaTValue()/2702.0/949.0/0.0012 + Tw_[i];
        //     }
        //     Info << "max Tw " << tab << db().time().timeName() << tab << gMax(Tw_) << endl;
        }
        else{
            UIndirectList<scalar>(Tw_, pp.faceCells()) =
            		T_.boundaryField()[patchi];
            // set char fraction internal field equal to boundary values
            UIndirectList<scalar>(charFrac_, pp.faceCells()) =
            		charFrac_.boundaryField()[patchi];
            // set pyrolysis radiative heat flux internal field to boundary values
            UIndirectList<scalar>(qemPyro_, pp.faceCells()) =
                        qemPyro_.boundaryField()[patchi];
        }
    }
    Tw_.correctBoundaryConditions();

    // Update heat transfer to wall (used in film/pyrolysis coupling)
    // heat flow out of film is positive
    qFilmToWall_ = htcw_->h()*(T_ - Tw_);

    // Info << "Tp " << T_[10] << endl;
    // Info << "Twp " << Tw_[10] << endl;
    // Info << "htc " << htcw_->h() << endl;

    forAll(qFilmToWall_,i){
        if(delta_[i]<1e-8){
            qFilmToWall_[i]=0.0;
        }
    }

    //scalarField qw(htcw_->h()*(T_-Tw_));
    //scalarField qs(htcs_->h()*(T_-TPrimary_));
    //DEBUG(qw[2]);
    //DEBUG(qs[2]);

    qFilmToWall_.correctBoundaryConditions();

	transmissivity_ = radiation_->transmissivity();
	transmissivity_.correctBoundaryConditions();
    // Update heat transfer from gas phase (used in diagnostics)
    // heat flow out of film is positive
    qGasToFilm_ = htcs_->h()*(T_ - TPrimary_);
    forAll(qGasToFilm_,i){
        if(delta_[i]<1e-8){
            qGasToFilm_[i]=0.0;
        }
    }
    qGasToFilm_.correctBoundaryConditions();

    // Update film surface temperature
    // TODO: Ts estimation should go here

    const scalarField d2k(delta_/2.0/kappa_);
    const scalarField Shs(radiation_->ShsConst());
    const scalarField hInf(htcs_->h());

    // Ts_ = 
    //     (T_+d2k*(hInf*TPrimary_+Shs))
    //     /(1+d2k*hInf);

    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hsSpPrimaryBf =
        hsSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hsSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();

    // Apply enthalpy source as difference between incoming and actual states
    hsSp_ -= rhoSp_*hs_;
}


void thermoSingleLayer::correctAlpha()
{
    /*DEBUG(hydrophilic_)*/
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos0(delta_ - dimensionedScalar("deltaWet", dimLength, deltaWet_));
    }
}


void thermoSingleLayer::updateSubmodels()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update vaporization
    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    //for diagnostics
    primaryMassTrans_.correctBoundaryConditions();
    primaryEnergyTrans_.correctBoundaryConditions();
    charFrac_.correctBoundaryConditions();
    qemPyro_.correctBoundaryConditions();

    // Update massAbsorption
    massAbsorption_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        massAbs_,
        energyAbs_
    );

    //for diagnostics
    massAbs_.correctBoundaryConditions();
    energyAbs_.correctBoundaryConditions();

    // Update radiation
    radiation_->correct();

    // Update kinematic sub-models
    kinematicSingleLayer::updateSubmodels();

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryEnergyTrans_);

    // Update source fields (phase change)
    hsSp_ += primaryEnergyTrans_/magSf()/time().deltaT();
    rhoSp_ += primaryMassTrans_/magSf()/time().deltaT();

    // Update source fields (mass absorption)
    rhoSp_ += massAbs_/magSf()/time().deltaT();

    // Vapour recoil pressure (can become unstable for wild oscillations in vaporization rate)
    // pSp_ -= sqr(primaryMassTrans_/magSf()/time_.deltaT())/2.0/rhoPrimary_;
    // Info << "vaporRecoilPressure " << gMin(pSp_) << " " << gAverage(pSp_) << " " << gMax(pSp_) << nl;
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& hs) const
{
    dimensionedScalar Tstd("Tstd", dimTemperature, 298.15);

    Info << "htcs_->h(): " << tab << db().time().timeName() << tab << max(htcs_->h()) << endl;
    Info << "htcw_->h(): " << tab << db().time().timeName() << tab << max(htcw_->h()) << endl;

    scalar maxTw = max(Tw_.primitiveField());
    reduce(maxTw,maxOp<scalar>());
    Info << "max(Tw_): " << tab << db().time().timeName() << tab << maxTw << endl;

    scalar maxTPrimary = max(TPrimary_.primitiveField());
    reduce(maxTPrimary,maxOp<scalar>());
    Info << "max(TPrimary_): " << tab << db().time().timeName() << tab << maxTPrimary << endl;

    return
       (
        - fvm::Sp(htcs_->h()/Cp_, hs) - htcs_->h()*(Tstd - TPrimary_)
        - fvm::Sp(htcw_->h()/Cp_, hs) - htcw_->h()*(Tstd - Tw_)
        );
}


void thermoSingleLayer::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    updateSurfaceTemperatures();

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
      // is vaporization energy accounted for twice?  hsSp_ and rhoSp_*hs_?
      - hsSp_
      - rhoSp_*hs_
      + q(hs_)
      //<luwi--option to include radiative emission from pyrolysis region>
      + radiation_->Shs(qemPyro_)
      //+ radiation_->Shs()
      //</luwi>
      // - fvm::SuSp(rhoSp_, hs_)
    );

    correctThermoFields();

    // Evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, regionType, false),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    Cp_
    (
        IOobject
        (
            "Cp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("Cp", dimEnergy/dimMass/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "kappa",
            dimEnergy/dimTime/dimLength/dimTemperature,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    T_
    (
        IOobject
        (
            "Tf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    hs_
    (
        IOobject
        (
            "hf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero),
        hsBoundaryTypes()
    ),
    qGasToFilm_
    (
     IOobject
     (
      "qGasToFilm",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qFilmToWall_
    (
     IOobject
     (
      "qFilmToWall",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),

    massAbs_
    (
        IOobject
        (
            "massAbsorption",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    energyAbs_
    (
        IOobject
        (
            "energyMassAbsorption",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    primaryEnergyTrans_
    (
        IOobject
        (
            "primaryEnergyTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    charFrac_
    (
        IOobject
        (
            "charFrac",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    //luwi
    qemPyro_
    (
        IOobject
        (
            "qemPyro",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    // luwi
    transmissivity_
    (
        IOobject
        (
            "filmTransmissivity",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    deltaWet_(coeffs_.get<scalar>("deltaWet")),
    hydrophilic_(coeffs_.get<bool>("hydrophilic")),
    hydrophilicDryScale_(0.0),
    hydrophilicWetScale_(0.0),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimArea/dimTime, Zero),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // Must have same name as hSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    YPrimary_(),

    viscosity_(filmViscosityModel::New(*this, coeffs(), mu_)),
    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    massAbsorption_(massAbsorptionModel::New(*this, coeffs())), // kvm
    radiation_(filmRadiationModel::New(*this, coeffs())),
    Tmin_(-VGREAT),
    Tmax_(VGREAT)
{
    if (coeffs().readIfPresent("Tmin", Tmin_))
    {
        Info<< "    limiting minimum temperature to " << Tmin_ << endl;
    }

    if (coeffs().readIfPresent("Tmax", Tmax_))
    {
        Info<< "    limiting maximum temperature to " << Tmax_ << endl;
    }

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, Zero),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    /*DEBUG(hydrophilic_)*/
    if (hydrophilic_)
    {
        coeffs_.readEntry("hydrophilicDryScale", hydrophilicDryScale_);
        coeffs_.readEntry("hydrophilicWetScale", hydrophilicWetScale_);
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctAlpha();

        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_);

        deltaRho_ == delta_*rho_;

        surfaceScalarField phi0
        (
            IOobject
            (
                "phi",
                time().timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            fvc::flux(deltaRho_*U_)
        );

        phi_ == phi0;

        // Evaluate viscosity from user-model
        viscosity_->correct(pPrimary_, T_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        Info<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void thermoSingleLayer::preEvolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

//    correctHsForMappedT();

    kinematicSingleLayer::preEvolveRegion();

    // Update phase change
    primaryMassTrans_ == dimensionedScalar("zero", dimMass, 0.0);
    primaryEnergyTrans_ == dimensionedScalar("zero", dimEnergy, 0.0);
    // Update mass absorption
    massAbs_ == dimensionedScalar("zero", dimMass, 0.0); // kvm
    energyAbs_ == dimensionedScalar("zero", dimEnergy, 0.0); // kvm
}


void thermoSingleLayer::evolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Update film coverage indicator
    correctAlpha();


    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update temperature using latest hs_
    T_ == T(hs_); // kvm
    //kvm-debug Info << "thermoSingleLayer::T_ " << T_ << endl;

#include "diagnostics.H"

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update film wall and surface velocities
    updateSurfaceVelocities(); // kvm

    // Update film wall and surface temperatures
    // updateSurfaceTemperatures();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const volScalarField& thermoSingleLayer::Cp() const
{
    return Cp_;
}


const volScalarField& thermoSingleLayer::kappa() const
{
    return kappa_;
}


const volScalarField& thermoSingleLayer::T() const
{
    return T_;
}


const volScalarField& thermoSingleLayer::Ts() const
{
    return Ts_;
}


const volScalarField& thermoSingleLayer::Tw() const
{
    return Tw_;
}


const volScalarField& thermoSingleLayer::hs() const
{
    return hs_;
}


tmp<volScalarField> thermoSingleLayer::massAbs() const // kvm
{
    return massAbs_; // kvm
}


tmp<volScalarField> thermoSingleLayer::primaryMassTrans() const
{
    return primaryMassTrans_;
}


void thermoSingleLayer::info()
{
    kinematicSingleLayer::info();

    const scalarField& Tinternal = T_;

    Info<< indent << "min/mean/max(T)    = "
        << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", "
        << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
    massAbsorption_->info(Info); // kvm
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    scalarField& Srho = tSrho.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchi = intCoupledPatchIDs()[i];

        scalarField patchMass =
            primaryMassTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchMass);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho(" + Foam::name(i) + ")",
                time_.timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );
/*
    phase change energy fed back into the film...

    scalarField& Sh = tSh.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }
*/
    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace surfaceFilmModels

// ************************************************************************* //
