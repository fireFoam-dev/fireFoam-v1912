/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "physicoChemicalConstants.H"
#include "fvDOM.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ThermoCloud<CloudType>&
Foam::ThermoCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ThermoCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ThermoCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::SLGThermo& Foam::ThermoCloud<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ThermoCloud<CloudType>::T() const
{
    return T_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ThermoCloud<CloudType>::p() const
{
    return p_;
}


template<class CloudType>
inline const Foam::HeatTransferModel<Foam::ThermoCloud<CloudType>>&
Foam::ThermoCloud<CloudType>::heatTransfer() const
{
    return *heatTransferModel_;
}


template<class CloudType>
inline const Foam::integrationScheme&
Foam::ThermoCloud<CloudType>::TIntegrator() const
{
    return *TIntegrator_;
}


template<class CloudType>
inline bool Foam::ThermoCloud<CloudType>::radiation() const
{
    return radiation_;
}

//ankur
template<class CloudType>
inline bool Foam::ThermoCloud<CloudType>::coupledRadiation() const
{
    return coupledRadiation_;
}

//ankur
template<class CloudType>
inline Foam::label Foam::ThermoCloud<CloudType>::nBands() const
{
    return nBands_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaP(const label bandI) // ankur
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaP_[bandI]; // ankur
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaP(const label bandI) const // ankur
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaP_[bandI]; // ankur
}

// ankur
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPSc(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPSc_[bandI];
}


// ankur
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPSc(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPSc_[bandI];
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radT4(const label bandI) // ankur
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radT4_[bandI]; // ankur
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radT4(const label bandI) const // ankur
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radT4_[bandI]; // ankur
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPT4(const label bandI) // ankur
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPT4_[bandI]; // ankur
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPT4(const label bandI) const // ankur
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPT4_[bandI]; // ankur
}

// ankur
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPScAsy(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPScAsy_[bandI];
}


// ankur
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::radAreaPScAsy(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPScAsy_[bandI];
}

// ankur
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::CScat(const label bandI)
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return CScat_[bandI];
}


// ankur
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::CScat(const label bandI) const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return CScat_[bandI];
}

// ankur
template<class CloudType>
inline const Foam::List<Foam::scalar>& Foam::ThermoCloud<CloudType>::diaVal() const
{
  return diaVal_;
}

// ankur
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsTrans()
{
    return hsTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsTrans() const
{
    return hsTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsCoeff()
{
    return hsCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<CloudType>::hsCoeff() const
{
    return hsCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermoCloud<CloudType>::Sh(volScalarField& hs) const
{
    if (debug)
    {
        Info<< "hsTrans min/max = " << min(hsTrans()).value() << ", "
            << max(hsTrans()).value() << nl
            << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
            << max(hsCoeff()).value() << endl;
    }

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("h"))
        {
            const volScalarField Cp(thermo_.thermo().Cp());
            const volScalarField::Internal
                Vdt(this->mesh().V()*this->db().time().deltaT());

            return
                hsTrans()/Vdt
              - fvm::SuSp(hsCoeff()/(Cp*Vdt), hs)
              + hsCoeff()/(Cp*Vdt)*hs;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -hsTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(hs, dimEnergy/dimTime);
}

//<fmglobal>
template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::Ep(const label bandI) const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    if (radiation_)
    {
        scalarField& Ep = tEp.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& sumAreaPT4 = radAreaPT4_[bandI].field(); // ankur

        //Ep = sumAreaPT4*epsilon*physicoChemical::sigma.value()/V/dt;
        // ankur, in this implementation, epsilon (i.e., emission efficiency) and division by V is included in radAreaP_ term
        Ep = sumAreaPT4*physicoChemical::sigma.value()/dt;
    }

    return tEp;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::ap(const label bandI) const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    if (radiation_)
    {
        scalarField& ap = tap.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& sumAreaP = radAreaP_[bandI].field(); // ankur

        //ap = sumAreaP*epsilon/V/dt;
        // ankur, in this implementation, epsilon (i.e., absorption efficiency) and division by V is included in radAreaP_ term
        ap = sumAreaP/dt;
    }

    return tap;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<CloudType>::sigmap(const label bandI) const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    if (radiation_)
    {
        scalarField& sigmap = tsigmap.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& sumAreaP = radAreaPSc_[bandI].field(); // ankur

        //sigmap *= sumAreaP*(1.0 - f)*(1.0 - epsilon)/V/dt;
        // ankur, in this implementation, scattering efficiency and division by V is included in radAreaP_ term
        sigmap = sumAreaP/dt;
    }

    return tsigmap;
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::phaseFunc(const label bandI, const label sour, const label dest) const
{
    tmp<volScalarField> tPhaseFunc
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:phaseFunc",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if ((radiation_) && (sour != dest))
    {
        scalarField& phaseF = tPhaseFunc.ref().primitiveFieldRef();
        const scalarField& cscat = CScat_[bandI].field();
        scalarField asyFacVal(radAreaPScAsy_[bandI].field()/(radAreaPSc_[bandI].field() + VSMALL));
        scalarField A_asy(2.*asyFacVal/(1. + asyFacVal*asyFacVal));
        scalarField B_asy((1. - asyFacVal*asyFacVal)/(1. + asyFacVal*asyFacVal));

        const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<radiation::radiationModel>("radiationProperties");
        const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);
        const Foam::vector dRayS = fvdomModel_.IRay(sour).d();
        const Foam::vector dRayD = fvdomModel_.IRay(dest).d();
        scalar sDotd = dRayS & dRayD;
        scalarField AsDotd(A_asy*sDotd);
#if 1
	phaseF =
	(
	    cscat * B_asy * (AsDotd*(AsDotd*(AsDotd*(AsDotd + 1.) + 1.) + 1.) + 1.)
	    / (4.*constant::mathematical::pi - fvdomModel_.IRay(sour).omega() )
	);
#else
        phaseF = cscat * B_asy * (1. + A_asy*sDotd + A_asy*A_asy*sDotd*sDotd + Foam::pow(A_asy,3.)*Foam::pow(sDotd,3.) + Foam::pow(A_asy,4.)*Foam::pow(sDotd,4.)) / (4.*constant::mathematical::pi - fvdomModel_.IRay(sour).omega() );
#endif
     }

    return tPhaseFunc;
}
//Henyey-Greenstein phase function, computed and normalized directly everytime
template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::phaseFuncHGDirect(const label bandI, const label sour, const label dest) const
{
    tmp<volScalarField> tPhaseFunc
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:phaseFuncHGDirect",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

//Compute Henyey Greenstein phase function
    scalarField& phaseF = tPhaseFunc.ref().primitiveFieldRef();
    scalarField asyFacVal(radAreaPScAsy_[bandI].field()/(radAreaPSc_[bandI].field() + VSMALL));
    scalarField asyFacVal2(asyFacVal*asyFacVal);
    const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<radiation::radiationModel>("radiationProperties");
    const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);
    const Foam::vector dRayS = fvdomModel_.IRay(sour).d();
    const Foam::vector dRayD = fvdomModel_.IRay(dest).d();
    scalar sDotd = dRayS & dRayD;
    scalar nRays = fvdomModel_.nRay();

    phaseF = (1. - asyFacVal2)/Foam::pow((1. + asyFacVal2 - 2.*asyFacVal*sDotd),1.5);
    // Then, normalize phase function
    scalarField sum=phaseF;
    for(label iRay=0; iRay<nRays; iRay++)
    {
	if( iRay == sour )
	    continue;
	scalar siDotd = fvdomModel_.IRay(iRay).d() & dRayD;
	sum += (1. - asyFacVal2)/Foam::pow((1. + asyFacVal2 - 2.*asyFacVal*siDotd),1.5);
    }
    phaseF /= sum;

    return tPhaseFunc;
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::ThermoCloud<CloudType>::asyFac(const label bandI) const
{
    tmp<volScalarField> tasyFac
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:asyFac",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    scalarField& asyFac = tasyFac.ref().primitiveFieldRef();

    asyFac = radAreaPScAsy_[bandI].field()/(radAreaPSc_[bandI].field() + VSMALL);

    return tasyFac;
}

template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::radCalc(const label celli, const scalar pDia, const scalar ap, const scalar np, const scalar T4, const scalar dt)
{
        const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<radiation::radiationModel>("radiationProperties");
        const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);

	scalar sigmaAbsG = 0.;

        List<scalar> QAbsorp(nBands_,0.);
        List<scalar> QScat(nBands_,0.);
        List<scalar> asyFactor(nBands_,0.);

        scalar npDen = np/this->mesh().V()[celli];

        getQAbsQSct(pDia,npDen,QAbsorp,QScat,asyFactor);

        scalar multFac = dt*npDen*ap;
        forAll(radAreaP_, bandI)
        {
            scalar QAbs = QAbsorp[bandI];
            scalar QSct = QScat[bandI];
            scalar asy  = asyFactor[bandI];
            radAreaP_[bandI].field()[celli] += multFac*QAbs;
            radAreaPSc_[bandI].field()[celli] += multFac*QSct;
            radT4_[bandI].field()[celli] += dt*npDen*T4;
            radAreaPScAsy_[bandI].field()[celli] += multFac*QSct*asy;

            sigmaAbsG += QAbs*fvdomModel_.GLambda(bandI)[celli];
        }

        return (sigmaAbsG *= ap);
}

template<class CloudType>
inline void Foam::ThermoCloud<CloudType>::getQAbsQSct(const scalar pDia, const scalar np, List<scalar>& QAbs, List<scalar>& QSct, List<scalar>& asyFac) const
{

    label begin = 0;
    label end = diaVal_.size()-1;
    label nDias = end + 1;

    List<scalar> qabs(numPropBands_,0.);
    List<scalar> qsct(numPropBands_,0.);
    List<scalar> asyF(numPropBands_,0.);

    if(nDias==1)
    {
        qabs = absEff_[begin];
        qsct = sctEff_[begin];
        asyF = asyFac_[begin];
    }
    else if(pDia <= diaVal_[begin])
    {
        qabs = (absEff_[begin+1]-absEff_[begin])/(diaVal_[begin+1]-diaVal_[begin])*(pDia - diaVal_[begin])+absEff_[begin];
        qsct = (sctEff_[begin+1]-sctEff_[begin])/(diaVal_[begin+1]-diaVal_[begin])*(pDia - diaVal_[begin])+sctEff_[begin];
        asyF = (asyFac_[begin+1]-asyFac_[begin])/(diaVal_[begin+1]-diaVal_[begin])*(pDia - diaVal_[begin])+asyFac_[begin];
    }
    else if(pDia >= diaVal_[end])
    {
        qabs = (absEff_[end]-absEff_[end-1])/(diaVal_[end]-diaVal_[end-1])*(pDia - diaVal_[end])+absEff_[end];
        qsct = (sctEff_[end]-sctEff_[end-1])/(diaVal_[end]-diaVal_[end-1])*(pDia - diaVal_[end])+sctEff_[end];
        asyF = (asyFac_[end]-asyFac_[end-1])/(diaVal_[end]-diaVal_[end-1])*(pDia - diaVal_[end])+asyFac_[end];
    }
    else
    {
        for(label i=0;i<diaVal_.size()-1;i++)
        {
            label i1 = i;
            label i2 = i+1;
            if(pDia > diaVal_[i1] && pDia <= diaVal_[i2])
            {
                scalar weight = (pDia-diaVal_[i1])/(diaVal_[i2]-diaVal_[i1]);
                qabs = (1.0-weight)*absEff_[i1]+weight*absEff_[i2];
                qsct = (1.0-weight)*sctEff_[i1]+weight*sctEff_[i2];
                asyF = (1.0-weight)*asyFac_[i1]+weight*asyFac_[i2];
            }
        }
    }

    if (nBands_ !=1)
    {
        QAbs = qabs;
        QSct = qsct;
        asyFac = asyF;
    }
    else
    {
        List<scalar> tau(numPropBands_,0.);
        scalar tauGray = 0.;
        scalar rSum = 0.;

        scalar asyNum = 0.;
        scalar asyDen = 0.;

        forAll(tau,pBandI)
        {
            tau[pBandI] = 1 - exp(-0.25*3.1415*pDia*pDia*np*(qabs[pBandI]+qsct[pBandI])*beamLen_);
            tauGray += energyFrac_[pBandI]*tau[pBandI];
            rSum += energyFrac_[pBandI]*tau[pBandI]*qabs[pBandI]/(qabs[pBandI]+qsct[pBandI]+VSMALL);
            asyNum += energyFrac_[pBandI]*qsct[pBandI]*asyF[pBandI];
            asyDen += energyFrac_[pBandI]*qsct[pBandI];
        }

        scalar kextGray = -log(1. - tauGray)/(beamLen_+VSMALL);
        QAbs[0] = kextGray*rSum/(tauGray+VSMALL);
        QSct[0] = kextGray - QAbs[0];
        scalar kOverQ = 0.25*3.1415*pDia*pDia*np;
        QAbs[0] = QAbs[0]/(kOverQ + VSMALL);
        QSct[0] = QSct[0]/(kOverQ + VSMALL);
        asyFac[0] = asyNum/(asyDen + VSMALL);
    }

    return;
}
//</fmglobal>

template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::Tmax() const
{
    scalar val = -GREAT;
    bool nonEmpty = false;

    for (const parcelType& p : *this)
    {
        val = max(val, p.T());
        nonEmpty = true;
    }

    if (returnReduce(nonEmpty, orOp<bool>()))
    {
        return returnReduce(val, maxOp<scalar>());
    }

    return 0;
}


template<class CloudType>
inline Foam::scalar Foam::ThermoCloud<CloudType>::Tmin() const
{
    scalar val = GREAT;
    bool nonEmpty = false;

    for (const parcelType& p : *this)
    {
        val = min(val, p.T());
        nonEmpty = true;
    }

    if (returnReduce(nonEmpty, orOp<bool>()))
    {
        return returnReduce(val, minOp<scalar>());
    }

    return 0;
}

// ankur
template<class CloudType>
void Foam::ThermoCloud<CloudType>::setCScatCoeffs()
{
     const Foam::radiation::radiationModel& radModel_ = this->mesh().template lookupObject<Foam::radiation::radiationModel>("radiationProperties");
     const Foam::radiation::fvDOM& fvdomModel_=refCast<const Foam::radiation::fvDOM>(radModel_);

     label nRays = fvdomModel_.nRay();

     CScatCoeffs_ = 0.;
     const Foam::vector dRay0 = fvdomModel_.IRay(0).d();
     scalar omega0 = fvdomModel_.IRay(0).omega();

     for (int i=1; i < nRays; i++)
     {
         scalar s0Dotsi = dRay0 & fvdomModel_.IRay(i).d();
         scalar omegaI  = fvdomModel_.IRay(i).omega();

         CScatCoeffs_[0] += omegaI*s0Dotsi;
         CScatCoeffs_[1] += omegaI*s0Dotsi*s0Dotsi;
         CScatCoeffs_[2] += omegaI*s0Dotsi*s0Dotsi*s0Dotsi;
         CScatCoeffs_[3] += omegaI*s0Dotsi*s0Dotsi*s0Dotsi*s0Dotsi;
     }

     CScatCoeffs_[0] *= 1./(4.*constant::mathematical::pi - omega0);
     CScatCoeffs_[1] *= 1./(4.*constant::mathematical::pi - omega0);
     CScatCoeffs_[2] *= 1./(4.*constant::mathematical::pi - omega0);
     CScatCoeffs_[3] *= 1./(4.*constant::mathematical::pi - omega0);
}

// ankur
template<class CloudType>
void Foam::ThermoCloud<CloudType>::setCScat()
{
    forAll(CScat_, bandI)
    {
        scalarField& cscat = CScat_[bandI].field();
        scalarField asyFacVal(radAreaPScAsy_[bandI].field()/(radAreaPSc_[bandI].field() + VSMALL));
        scalarField A_asy(2.*asyFacVal/(1. + asyFacVal*asyFacVal));
        scalarField B_asy((1. - asyFacVal*asyFacVal)/(1. + asyFacVal*asyFacVal));
        cscat =  1./( B_asy*(1. + A_asy*CScatCoeffs_[0] + A_asy*A_asy*CScatCoeffs_[1] + Foam::pow(A_asy,3.)*CScatCoeffs_[2] + Foam::pow(A_asy,4.)*CScatCoeffs_[3]) + VSMALL);
    }
}

//<fmglobal>
template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<CloudType>::sprayTemperature() const
{
    tmp<volScalarField> tsprayTemperature
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":sprayTemperature",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTemperature, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    tmp<volScalarField> tsprayMass
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":sprayMass",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& sprayTemperature = tsprayTemperature.ref();
    volScalarField& sprayMass = tsprayMass.ref();
    for (const parcelType& p : *this)
    {
        const label celli = p.cell();

        sprayTemperature[celli] += p.nParticle()*p.mass()*p.T();
        sprayMass[celli] += p.nParticle()*p.mass();
    }

    forAll(sprayMass,celli)
    {
        if(sprayMass[celli] <= 0)
        {
            sprayTemperature[celli] = 293.15;
        }
        else
        {
            sprayTemperature[celli] /= sprayMass[celli];
        }
    }
    //sprayTemperature.primitiveFieldRef() /= mesh_.V();
    sprayTemperature.correctBoundaryConditions();

    return tsprayTemperature;
}
//</fmglobal>

// ************************************************************************* //
