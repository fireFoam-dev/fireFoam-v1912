void thermoSingleLayerPw::computeContactAngle() //called from initialise
{
    if (debug)
    {
        tabAdd();
        INFO<< "thermoSingleLayerPw::computeContactAngel()" << endl;
    }

    static volScalarField randomField
        (
         IOobject
         (
          "randomField",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
        );
    static volScalarField randomField2
        (
         IOobject
         (
          "randomField2",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
        );
    static label first=1;
    if(first){
        static label seed=0;
        Random ranGen_(seed);//kvm
        forAll(randomField,index){
            //normal distribution
            ranGen_.randomise01(randomField[index]);
            //normal distribution
            ranGen_.randomise01(randomField2[index]);
            /*Box-Muller approach for generating a random number from a normal distribution*/
            /*=SQRT(-2*LN(RAND()))*SIN(2*PI()*RAND())*/
            randomField[index]=sqrt(-2.0*log(randomField[index]))*sin(2.0*constant::mathematical::pi*randomField2[index]);
        }
        first=0;
    }

    forAll(contactAngle_,index){
        // normal distribution
        contactAngle_[index]=(randomField[index]*contactAngleStdDev_+contactAngleMean_);
        contactAngle_[index]=max(contactAngle_[index],contactAngleMin_);
        contactAngle_[index]=min(contactAngle_[index],contactAngleMax_);
    }

    if (debug)
    {
        INFO<< "leaving thermoSingleLayerPw::computeContactAngel()" << endl;
        tabSubtract();
    }
    return;
}

void thermoSingleLayerPw::updateOmega()
{
    if (debug)
    {
        tabAdd();
        INFO<< "thermoSingleLayerPw::updateOmega()" << endl;
    }
    scalar wet=1.0;
    scalar dry=0.0;
    /*DEBUG(hydrophilic_)*/
    if(!hydrophilic_){
        omega_=pos(delta_-wetToggle_*criticalFilmThickness_)*wet;
    }
    else{
        forAll(omega_,index){
            /* If film thickness is greater than the minimum film thickness, then set face to wet */
            if(omega_[index]==dry&&delta_[index] > wetToggle_*criticalFilmThickness_.value()){
                omega_[index]=wet;
                timeWetted_[index]=time_.time().value();
            }
            /* If film thickness is much much less than the minimum film thickness, then set face to dry */
            else if(omega_[index]==wet&&delta_[index]<dryToggle_*criticalFilmThickness_.value()){
                omega_[index]=dry;
            }
        }
    }
    omega_.correctBoundaryConditions();

    if (debug)
    {
        INFO<< "leaving thermoSingleLayerPw::updateOmega()" << endl;
        tabSubtract();
    }
    return;
}

void thermoSingleLayerPw::updateAlpha()
{
    /*DEBUG("setting alpha_ = omega_");*/
    alpha_ = omega_;
    return;
}

void thermoSingleLayerPw::initialise()
{
    if (debug)
    {
        tabAdd();
        INFO<< "thermoSingleLayerPw::initialise()" << endl;
    }

    read();

    /*if contact angle read from file, then don't compute*/
    if(!(max(contactAngle_).value()>0)){
        contactAngleFromFile_=false;
        computeContactAngle();
    }
    else{
        contactAngleFromFile_=true;
    }

    if (debug)
    {
        INFO<< "leaving thermoSingleLayerPw::initialise()" << endl;
        tabSubtract();
    }
}

void thermoSingleLayerPw::updateContactLine()
{
    if (debug)
    {
        tabAdd();
        INFO<< "thermoSingleLayerPw::updateContactLine()" << endl;
    }

    volVectorField gradOmega(fvc::grad(omega_));

    /*eliminate unnecessary noise in contactLine calculation*/
    /*0.1 is an arbitrary choice.  Any small number will do.  (1-0)/dx, as long as dx is less than 10 */
    forAll(gradOmega,index){
        if(mag(gradOmega[index])<0.1){
            gradOmega[index]=vector::zero;
        }
    }

    contactLine_=(gradOmega/(mag(gradOmega)+dimensionedScalar("SMALL",pow(dimLength,-1),SMALL)));

    // Remove contact line in cells that are wet
    //contactLine_*=(1-omega_);

    // Remove any patch-normal components of contactLine
    contactLine_ -= nHat()*(nHat() & contactLine_);
    if (debug)
    {
        INFO<< "leaving thermoSingleLayerPw::updateContactLine()" << endl;
        tabSubtract();
    }
}

const volScalarField& thermoSingleLayerPw::omega() const
{
    return omega_;
}

// deprecated const volScalarField& thermoSingleLayerPw::Yv() const
// deprecated {
// deprecated     return Yv_;
// deprecated }

void thermoSingleLayerPw::updateMassAbsorption()
{
    if (debug)
    {
        tabAdd();
        INFO<< "thermoSingleLayerPw::updateMassAbsorption()" << endl;
    }
    const scalar T1=20.0+273.15; //K
    const scalar A1=0.0207; //kg/m2
    const scalar n1=0.456;
    const scalar T2=43.0+273.15; //K
    const scalar A2=0.0317; //kg/m2
    const scalar n2=0.465;
    const scalar dt=time_.deltaT().value();

    forAll(massAbsorption_,index){
        if(omega_[index]==1){
            /*linear interpolation for absorption pre-factor*/
            scalar Tloc=max(T()[index],T1);
            Tloc=min(Tloc,T2);
            scalar A = (Tloc-T1)/(T2-T1)*(A2-A1)+A1;
            scalar n = (Tloc-T1)/(T2-T1)*(n2-n1)+n1;
            const scalar t=max(pow(cummulativeMassAbsorption_[index]/A,1/n)*60.0,dt);
            const scalar to=max(SMALL,t-dt);
            massAbsorption_[index]=A*(pow(t/60.0,n)-pow((to)/60.0,n))*magSf()[index];
            cummulativeMassAbsorption_[index]+=massAbsorption_[index]/magSf()[index];
        }
        else{
            massAbsorption_[index]=0.0;
        }
    }
    if (debug)
    {
        INFO<< "leaving thermoSingleLayerPw::updateMassAbsorption()" << endl;
        tabSubtract();
    }
}


