//    //For diagnostics
//    if(diagnostics_){
//        sensibleEnthalpyChange_=delta_*rho_*(hs_-hs_.oldTime())/time().deltaT()+hs_*(deltaRho_-deltaRho_.oldTime())/time().deltaT();
//        sensibleEnthalpyChange_.correctBoundaryConditions();
//        //for post processing, need to divide by thickness of extruded mesh to get J/s
//        advectiveEnthalpy_=fvc::interpolate(hs_)*phi_;//
//        //advectiveEnthalpy_=delta_*rho_* (U_ & gTan()) * hs_;//phi_*hs_;//
//        //advectiveEnthalpy_.correctBoundaryConditions();
//    }
//
//    //diagnostics
//    if(diagnostics_){
//        dimensionedScalar Tb("Tboil", dimTemperature, 374.8);
//        dimensionedScalar Tinlet("Tinlet", dimTemperature, 298.15);
//        dimensionedScalar Patm("Patm", dimPressure, 101325.0);
//        dimensionedScalar Tstd("Tstd", dimTemperature, 298.15);
//        const liquid& liq = thermo_.liquids().properties()[liquidId_];
//        dimensionedScalar Cpb("Cpb", dimEnergy/dimMass/dimTemperature, liq.cp(Patm.value(), Tb.value()));
//        dimensionedScalar Cpi("Cpi", dimEnergy/dimMass/dimTemperature, liq.cp(Patm.value(), Tinlet.value()));
//        /*Info << "Cp " << Cpi << endl;*/
//
//        energyPhaseChangeDiagnostic_ = primaryEnergyTrans_/magSf()/time().deltaT() + primaryMassTrans_/magSf()/time().deltaT()*hs_ 
//        //energyPhaseChangeDiagnostic_ = primaryEnergyTrans_/magSf()/time().deltaT() + primaryMassTrans_/magSf()/time().deltaT()*hs_ //works with wallVaporization case
//        //energyPhaseChangeDiagnostic_ = primaryEnergyTrans_/magSf()/time().deltaT() //+ massPhaseChangeForPrimary_/magSf()/time().deltaT()*hs_
//            ;//- (Cpb*(Tb-Tstd)-Cpi*(Tinlet-Tstd))*massPhaseChangeForPrimary_/magSf()/time().deltaT(); //sensible enthalpy portion
//        energyPhaseChangeDiagnostic_.correctBoundaryConditions();
//    }

    if(diagnostics_){
        static dimensionedScalar qRad("qRad", dimEnergy, 0.0);
        qRad+=sum(radiation_->Shs()*magSf()*time_.deltaT());   
//        Info << "qRad " << time_.value() << " " << qRad.value() << endl;
        char buffer[256];
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qRad",time_.value(),qRad.value());
        Info << buffer;
        static dimensionedScalar qConvWall("qConvWall", dimEnergy, 0.0);
        dimensionedScalar small("SMALL", delta_.dimensions(), 0.0001);
        dimensionedScalar smallDeltaRho("SMALL", deltaRho_.dimensions(), 0.000000001);
        qConvWall-=sum(htcw_->h()*(T_ - Tw_)*magSf()*time_.deltaT());   
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qConvWall",time_.value(),qConvWall.value());
        Info << buffer;
        static dimensionedScalar qConvGas("qConvGas", dimEnergy, 0.0);
        qConvGas-=sum(deltaRho_*htcs_->h()*(T_ - TPrimary_)*magSf()*time_.deltaT())/sum(deltaRho_+smallDeltaRho);   
//        Info << "qConvGas " << time_.value() << " " << qConvGas.value() << endl;
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qConvGas",time_.value(),qConvGas.value());
        Info << buffer;
//        static dimensionedScalar qPhaseChange("qPhaseChange", dimEnergy, 0.0);
//        qPhaseChange-=sum(primaryEnergyTrans_);   
////        Info << "qPhaseChange " << time_.value() << " " << qPhaseChange.value() << endl;
//        sprintf(buffer,"%15s %10.5f % 10.5e\n","qPhaseChange",time_.value(),qPhaseChange.value());
//        Info << buffer;
        static dimensionedScalar qSensible("qSensible", dimEnergy, 0.0);
        //qSensible-=sum(deltaRho_*(hs_-hs_.oldTime())*magSf());   
        qSensible-=sum((delta_*rho_*hs_-delta_.oldTime()*rho_.oldTime()*hs_.oldTime())*magSf());   
//        Info << "qSensible " << time_.value() << " " << qSensible.value() << endl;
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qSensible",time_.value(),qSensible.value());
        Info << buffer;
        static dimensionedScalar qImp("qImp(rhoSp*hs)", dimEnergy, 0.0);
        qImp-=sum(rhoSp_*hs_.oldTime()*magSf()*time_.deltaT());   
//        Info << "qImp " << time_.value() << " " << qImp.value() << endl;
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qImp(rhoSp*hs)",time_.value(),qImp.value());
        Info << buffer;
        static dimensionedScalar qImpSens("qImpSens(hsSp_)", dimEnergy, 0.0);
        qImpSens-=sum(hsSp_*magSf()*time_.deltaT());   
//        Info << "qImpSens " << time_.value() << " " << qImpSens.value() << endl;
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qImpSens(hsSp_)",time_.value(),qImpSens.value());
        Info << buffer;
        static dimensionedScalar qTotal("qTotal", dimEnergy, 0.0);
        qTotal=qRad+qConvWall+qConvGas+qSensible+qImpSens+qImp;
//        Info << "qTotal " << time_.value() << " " << qConvWall.value()+qConvGas.value()+qPhaseChange.value()+qSensible.value()+qImpSens.value()+qImp.value() << endl;
        sprintf(buffer,"%15s %10.5f % 10.5e\n","qTotal",time_.value(),qTotal.value());
        Info << buffer;
//        Info << "maxDeltaf " << max(delta_).value() << endl;
//        static dimensionedScalar qAdv("qAdv", dimEnergy, 0.0);
//        qAdv-=sum(fvc::surfaceSum(phi_*fvc::interpolate(hs_)));   
//        Info << "qAdv " << time_.value() << " " << qAdv.value() << endl;
    }
    
