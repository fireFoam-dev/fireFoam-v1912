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

#include "thermoSingleLayerFmPanel.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedWallPolyPatch.H"
#include "constants.H" 

// Sub-models
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"
#include "stdio.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace regionModels
    {
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(thermoSingleLayerFmPanel, 0);
        addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerFmPanel, mesh);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

/*called from surfaceFilmModel::evolve()*/
bool thermoSingleLayerFmPanel::read()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerFmPanel::read()" << endl;
    }
    // no additional properties to read
    if(thermoSingleLayerPw::read())
    {
        augmentedRadiation_ = coeffs_.lookupOrDefault<Switch>("augmentedRadiation",false);
        if(augmentedRadiation_){
            const dictionary& subdict=coeffs_.subDict("augmentedRadiationCoeffs");
            subdict.lookup("qRadConstant") >> qRadConstant_;
            subdict.lookup("qRadXMax") >> qRadXMax_; 
            subdict.lookup("qRadXMin") >> qRadXMin_; 
            subdict.lookup("qRadYMax") >> qRadYMax_; 
            subdict.lookup("qRadYMin") >> qRadYMin_; 
            subdict.lookup("qRadBegin") >> qRadBegin_; 
            subdict.lookup("qRadEnd") >> qRadEnd_; 
            subdict.lookup("qRadEmissivity") >> qRadEmissivity_; 
            subdict.lookup("qRadAbsorptivity") >> qRadAbsorptivity_; 
        }

        solveLumpedCapacitance_ = coeffs_.lookupOrDefault<Switch>("solveLumpedCapacitance",false);
        if(solveLumpedCapacitance_){
            const dictionary& subdict=coeffs_.subDict("solveLumpedCapacitanceCoeffs");
            subdict.lookup("thickness") >> (lc.thickness); //thickness of aluminum panel
            subdict.lookup("density") >> lc.density;
            subdict.lookup("cp1") >> lc.cp1;
            subdict.lookup("cp2") >> lc.cp2;
            subdict.lookup("k1") >> lc.k1;
            subdict.lookup("k2") >> lc.k2;
            subdict.lookup("T1") >> lc.T1;
            subdict.lookup("T2") >> lc.T2;
            subdict.lookup("Tinit") >> lc.Tinit;
            subdict.lookup("hConvBack") >> lc.hConvBack;
            subdict.lookup("hConvFront") >> lc.hConvFront;
            subdict.lookup("TsurrBack") >> lc.TsurrBack;
            subdict.lookup("TsurrFront") >> lc.TsurrFront;
            subdict.lookup("emissivityFront") >> lc.emissivityFront;
            subdict.lookup("emissivityBack") >> lc.emissivityBack;
            subdict.lookup("convectiveScaling") >> lc.convectiveScaling;
        }
        
        xiangyang_ = coeffs_.lookupOrDefault<Switch>("XiangYang",false);
        Info << "looking up perfectlyWettedInlet\n";
        perfectlyWettedInlet_ = coeffs_.lookupOrDefault<Switch>("perfectlyWettedInlet",false);
        if(perfectlyWettedInlet_){
            const dictionary& subdict=coeffs_.subDict("perfectlyWettedInletCoeffs");
//            perfectlyWettedInletName_ = subdict.lookup("inletName");
            const wordList wordListTmp
            (
            	subdict.lookup("inletNames")
            );
            perfectlyWettedInletNames_=wordListTmp;
            subdict.lookup("offsetDistance") >> perfectlyWettedInletDistance_ ;
        }
        if (debug)
        {
            Info<<fmtab.c_str()<< "leaving thermoSingleLayerFmPanel::read()" << endl;
            tabSubtract();
        }
        return true;
    }
    else
    {
        return false;
    }
}

void thermoSingleLayerFmPanel::initialise()
{
    if (debug)
    {
        Pout<< "thermoSingleLayerFmPanel::initialise()" << endl;
    }

    read();

}

void thermoSingleLayerFmPanel::updateSubmodels()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerPw::updateSubmodels()" << endl;
    }

    thermoSingleLayerPw::updateSubmodels();

    /*partially wetted treatment*/
    if(perfectlyWettedInlet_){
        zeroContactAngleInlet();
    }
    thermoSingleLayerPw::updateContactLine();

    if(augmentedRadiation_){
        updateQRad();
    }

    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerPw::updateSubmodels()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayerFmPanel::updateSurfaceTemperatures()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerFmPanel::updateSurfaceTemperatures2()" << endl;
    }

    if(!solveLumpedCapacitance_){
        thermoSingleLayerPw::updateSurfaceTemperatures();
    }
    else{
        // Push boundary film temperature values into internal field
        static label first=1;
        if(first){
            for (label i=0; i<intCoupledPatchIDs_.size(); i++)
            {
                //label patchI = intCoupledPatchIDs_[i];
                //const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
//This doesn't seem to be very good for restarting.  It resets the wall temperature to the gas-phase boundary temperature
//                UIndirectList<scalar>(Tw_, pp.faceCells()) =
//                    T_.boundaryField()[patchI];
            }
            first=0;
        }

        /*lumped-capacitance bc model for aluminum plate*/
        if(solveLumpedCapacitance_){
            scalar filmDeltaWet=0.00002; //m
            scalar filmDeltaDry=0.0000; //m

            forAll(wettedFraction_,i){
                wettedFraction_[i]=(delta_[i]-filmDeltaDry)/(filmDeltaWet-filmDeltaDry);
            }
			wettedFraction_.max(0.0);
			wettedFraction_.min(1.0);

            qDotFilm_ = wettedFraction_*htcw_->h()*(T_ - Tw_);
            

            scalar g=9.8;
            scalar Tw=Tw_.average().value();
            scalar Tamb=294.0;
            scalar Tf=0.5*(Tw+Tamb);
            scalar beta=1.0/Tf;
            scalar L=0.9;
            scalar alpha=38.3e-6;
            scalar nu=26.4e-6;

            scalar RaL = g*beta*(Tw-Tamb)*pow(L,3)/(nu*alpha);

            scalar Pr=0.69;
            scalar NuL=pow(
                    0.825+(0.387)*pow(RaL,1./6.)/
                        pow(
                            1.0+pow(
                                0.492/Pr
                            ,9./16.)
                       ,8./27.)
                    ,2.);
            
            scalar k=33.8e-3;
            dimensionedScalar hbar("hbar",dimEnergy/dimArea/dimTime/dimTemperature,lc.convectiveScaling.value()*NuL*k/L);

            qDotBackDry_ = hbar*(lc.Tinit-Tw_);
            Info << "hbar = " << hbar << endl;
            Info << "RaL = " << RaL << endl;
            qDotBackDry_ = hbar*(lc.TsurrBack-Tw_);
            qDotFrontDry_ = (1.0-wettedFraction_)*hbar*(TPrimary_-Tw_);
            //qDotBackDry_ = lc.hConvBack*(lc.Tinit-Tw_);
            //qDotFrontDry_ = (1.0-wettedFraction_)*lc.hConvFront*(TPrimary_-Tw_);
            qRadAugmentedDry_ = (1.0-wettedFraction_)*qRadAugmented_;
            const dimensionedScalar dt = time().deltaT();

            #include "solveSolid.H"
        }

        qFilmToWall_ = htcw_->h()*(T_ - Tw_);
        qFilmToWall_.correctBoundaryConditions();
//        Info << "max htcw " << max(htcw_->h()) << endl;
//        Info << "max qFilmToWall_ " << max(qFilmToWall_) << endl;
//        Info << "min qFilmToWall_ " << min(qFilmToWall_) << endl;
        //Info << qFilmToWall_<<endl;
        
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
        Ts_ = T_;
        Ts_.correctBoundaryConditions();
    }
    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerFmPanel::updateSurfaceTemperatures2()" << endl;
        tabSubtract();
    }
}

tmp<fvScalarMatrix> thermoSingleLayerFmPanel::q
(
 volScalarField& hs
 ) const
{
    dimensionedScalar Tstd("Tstd", dimTemperature, 298.15);
    return
        (
         - fvm::Sp(htcs_->h()/Cp_, hs) - htcs_->h()*(Tstd - TPrimary_)
         - fvm::Sp(htcw_->h()/Cp_, hs) - htcw_->h()*(Tstd - Tw_)
         /* When the film becomes dry then qRadAugmented_ is to strong unless
          * heat loss (conduction through solid, surface emission, etc) 
          * is accounted for */
      //+ omega_*qRadAugmented_
      //TODO:  I am changing qRadAugmented_ here, but in standardPhaseChange qRadAugmented_ is unmodified
      + qRadAbsorptivity_*qRadAugmented_*wettedFraction_
      //this term can case temp to be wildly unstable!!!     - qRadEmissivity_*sigmaSB_*pow(T_,4) //kvm, need to lineraize this term and use hs/Cp_ instead of T_
    );
}

void thermoSingleLayerFmPanel::zeroContactAngleInlet()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerFmPanel::zeroContactAngleInlet()" << endl;
    }
    /*force contact angle to zero near inlet*/
    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();
    const volVectorField& cellCentres = regionMesh().C();
    static bool first=true;
    if(!contactAngleFromFile_){
        if(first){
            forAll(bm,patchI){

                for(label iName=0;iName<perfectlyWettedInletNames_.size();iName++)
                {
                	word perfectlyWettedInletName=perfectlyWettedInletNames_[iName];
					if(bm[patchI].name()== perfectlyWettedInletName){
						const pointField& bCenters = bm[patchI].faceCentres();
						forAll(cellCentres,i){
							scalar distance=1e10;
							forAll(bCenters,j){
								scalar tmpDistance=sqrt(
								   +pow(bCenters[j][0]-cellCentres[i][0],2)
								   +pow(bCenters[j][1]-cellCentres[i][1],2)
								   +pow(bCenters[j][2]-cellCentres[i][2],2)
								   );
								distance=(tmpDistance<distance)?(tmpDistance):(distance);
							}
							scalar distance1=perfectlyWettedInletDistance_.value();
							scalar distance2=2.0*distance1;
							if(distance<distance1){
								contactAngle_[i]=0.0;
							}
							else if(distance<distance2){
								contactAngle_[i]*=(distance-distance1)/(distance2-distance1);
							}
						}
					}
                }
            }
            first=false;
        }
    }

    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerFmPanel::zeroContactAngleInlet()" << endl;
        tabSubtract();
    }
}

void thermoSingleLayerFmPanel::updateQRad()
{
    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerFmPanel::updateQRad()" << endl;
    }

    const volVectorField& cellCentres = regionMesh().C();
    forAll(qRadAugmented_,index){
        if(time().time().value() >= qRadBegin_.value() && time().time().value() <= qRadEnd_.value()){
            if(cellCentres[index][1]>qRadYMin_.value()&&cellCentres[index][1]<qRadYMax_.value()){
                if(cellCentres[index][0]>qRadXMin_.value()&&cellCentres[index][0]<qRadXMax_.value()){
                    //assume uniform heat flux for now
                    //qRadAugmented_[index]=qRadConstant_.value()-5.67e-8*pow((T_[index]-300.0),4);
                    qRadAugmented_[index]=qRadConstant_.value();
                }
            }
        }
        else{
            qRadAugmented_[index]=0.0;
        }

    }
    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerFmPanel::updateQRad()" << endl;
        tabSubtract();
    }
}

void thermoSingleLayerFmPanel::wettedAreaInfo() const
{
    scalar wetArea=0.0;
    scalar dryArea=0.0;
    scalar totalArea=0.0;
    forAll(omega_,index){
        totalArea+=magSf()[index];
        if(omega_[index]==1){
            wetArea+=magSf()[index];
        }
        else if(omega_[index]==0){
            dryArea+=magSf()[index];
        }
    }

    const scalar totalWetArea=returnReduce(wetArea, sumOp<scalar>());
    const scalar totalDryArea=returnReduce(dryArea, sumOp<scalar>());
    Info << indent << "wettedAreaFraction  = " << totalWetArea/(totalDryArea+totalWetArea+ROOTVSMALL) << nl;
    Info << indent << "total area          = " << (totalDryArea+totalWetArea) << " m2"<< nl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerFmPanel::thermoSingleLayerFmPanel
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    thermoSingleLayerPw(modelType, mesh, g, regionType),
    augmentedRadiation_(false),
    xiangyang_(false),
    perfectlyWettedInlet_(false),
    perfectlyWettedInletNames_(
//        	subdict.lookup("inletNames")
	),
    perfectlyWettedInletDistance_(0.0),
    qRadAugmented_
    (
     IOobject
     (
      "qRadAugmented",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qRadAugmentedDry_
    (
     IOobject
     (
      "qRadAugmentedDry",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qDotFilm_
    (
     IOobject
     (
      "qDotFilm",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qDotBackDry_
    (
     IOobject
     (
      "qDotBack",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qDotFrontDry_
    (
     IOobject
     (
      "qDotFrontDry",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    wettedFraction_
    (
     IOobject
     (
      "wettedFraction",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimless, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
//    qFilmToWall_
//    (
//     IOobject
//     (
//      "qWall",
//      time().timeName(),
//      regionMesh(),
//      IOobject::NO_READ,
//      IOobject::AUTO_WRITE
//     ),
//     regionMesh(),
//     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0)
//    ),
    qRadConstant_(0.0),
    qRadXMax_(0.0),
    qRadXMin_(0.0),
    qRadYMax_(0.0),
    qRadYMin_(0.0),
    qRadBegin_(0.0),
    qRadEnd_(0.0),
    qRadEmissivity_(1.0),
    qRadAbsorptivity_(1.0),
    solveLumpedCapacitance_(false),
    lc()
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerFmPanel::~thermoSingleLayerFmPanel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayerFmPanel::postEvolveRegion()
{

    if (debug)
    {
        tabAdd();
        Info<<fmtab.c_str()<< "thermoSingleLayerFmPanel::postEvolveRegion()" << endl;
    }

    if (time().outputTime())
    {
        //massPhaseChangeForPrimary_.write();
        htcw_->h()->write();
        htcs_->h()->write();
        Tw_.write();
        qRadAugmented_.write();
        // scalarField& qWall = qWall().ref();
        qWall()().write(); // kvm, why do we need two sets of paranthesis?
#include "manualDecomposition.H"
    }
    if(diagnostics_){
        static dimensionedScalar qConvGas("qConvGas", dimEnergy, 0.0);
        qConvGas-=sum(omega_*htcs_->h()*(T_ - TPrimary_)*magSf()*time_.deltaT());   
        INFO << "qConvGas " << time_.value() << " " << qConvGas.value() << endl;
        static dimensionedScalar qPhaseChange("qPhaseChange", dimEnergy, 0.0);
        qPhaseChange-=sum(primaryEnergyTrans_);   
        INFO << "qPhaseChange " << time_.value() << " " << qPhaseChange.value() << endl;
        static dimensionedScalar qSensible("qSensible", dimEnergy, 0.0);
        //qSensible-=sum(delta_*rho_*(hs_-hs_.oldTime())*magSf());   
        qSensible-=sum((delta_*rho_*hs_-delta_.oldTime()*rho_.oldTime()*hs_.oldTime())*magSf());   
        INFO << "qSensible " << time_.value() << " " << qSensible.value() << endl;
        static dimensionedScalar qImp("qImp", dimEnergy, 0.0);
        qImp-=sum(rhoSp_*hs_*magSf()*time_.deltaT());   
        INFO << "qImp " << time_.value() << " " << qImp.value() << endl;
        static dimensionedScalar qImpSens("qImpSens", dimEnergy, 0.0);
        qImpSens-=sum(hsSp_*magSf()*time_.deltaT());   
        INFO << "qImpSens " << time_.value() << " " << qImpSens.value() << endl;
        INFO << "qTotal " << time_.value() << " " << qConvGas.value()+qPhaseChange.value()+qSensible.value()+qImpSens.value()+qImp.value() << endl;
//        static dimensionedScalar qAdv("qAdv", dimEnergy, 0.0);
//        qAdv-=sum(fvc::surfaceSum(phi_*fvc::interpolate(hs_)));   
//        INFO << "qAdv " << time_.value() << " " << qAdv.value() << endl;
    }
    if (xiangyang_){
        integrateSplashMass();
    }
    if (debug)
    {
        Info<<fmtab.c_str()<< "leaving thermoSingleLayerFmPanel::postEvolveRegion()" << endl;
        tabSubtract();
    }
}

tmp<DimensionedField<scalar, volMesh>>
thermoSingleLayerFmPanel::qWall() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
     qFilmToWall_
    );
}

tmp<DimensionedField<scalar, volMesh> >
thermoSingleLayerFmPanel::qRad() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        qRadAugmented_
    );
}

void thermoSingleLayerFmPanel::integrateSplashMass(){
    static const vector center(0,0,0);
    static scalarList binRadius(6,0.0);
    binRadius[0]=0.022;
    binRadius[1]=0.040;
    binRadius[2]=0.060;
    binRadius[3]=0.089;
    binRadius[4]=0.134;
    binRadius[5]=0.190;
    //const scalarList binRadius(.0,.018,.040,.060,.089,.134,.190);
    static const volVectorField& cellCentres = regionMesh().C();
    static scalarList cummulatedMass(binRadius.size(),0.0);

    forAll(delta_,i){
        scalar x=cellCentres[i][0];
        scalar y=cellCentres[i][1];
        scalar radius=sqrt(pow(x-center.x(),2)+pow(y-center.y(),2));
        for(label j=0;j<binRadius.size();j++){
            if(j==0){
                if(radius<binRadius[j]){
                    cummulatedMass[j]-=rhoSp_[i]*magSf()[i]*time().deltaT().value();
                }
            }
            else{
                if(binRadius[j-1]<=radius&&radius<binRadius[j]){
                    cummulatedMass[j]-=rhoSp_[i]*magSf()[i]*time().deltaT().value();
                }
            }

        }
    }
    static char buffer[256];
    sprintf(buffer,"SplashedMass %10.5f ",time().value());
    Info << buffer;
    forAll(cummulatedMass,j){
        sprintf(buffer," % 10.5g ",cummulatedMass[j]);
        Info << buffer;
    }
    sprintf(buffer," % 10.5g ",sum(cummulatedMass));
    Info << buffer;
    Info << endl;

}

void thermoSingleLayerFmPanel::info()
{
    thermoSingleLayerPw::info();

    wettedAreaInfo();

    if(solveLumpedCapacitance_){
        Info<< indent << "min/max(Tw)         = " << min(Tw_).value() << ", " << max(Tw_).value() << nl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam


// ************************************************************************* //
