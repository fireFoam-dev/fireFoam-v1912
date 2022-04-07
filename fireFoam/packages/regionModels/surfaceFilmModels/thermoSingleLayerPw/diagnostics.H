    if(diagnostics_){
        static dimensionedScalar qConvGas("qConvGas", dimEnergy, 0.0);
        qConvGas-=sum(htcs_->h()*(T_ - TPrimary_)*magSf()*time_.deltaT());   
        INFO << "qConvGas " << time_.value() << " " << qConvGas.value() << endl;
        static dimensionedScalar qPhaseChange("qPhaseChange", dimEnergy, 0.0);
        qPhaseChange-=sum(primaryEnergyPCTrans_);   
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
        INFO << "qTotal " << time_.value() << " " << qConvGas.value()+qPhaseChange.value()+qSensible.value()+qImpSens.value()+qImp.value() << endl;
        INFO << "maxDeltaf " << max(delta_).value() << endl;
//        static dimensionedScalar qAdv("qAdv", dimEnergy, 0.0);
//        qAdv-=sum(fvc::surfaceSum(phi_*fvc::interpolate(hs_)));   
//        Info << "qAdv " << time_.value() << " " << qAdv.value() << endl;
    }
