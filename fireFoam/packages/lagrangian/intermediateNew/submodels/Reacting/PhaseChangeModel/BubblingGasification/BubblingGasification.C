/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "BubblingGasification.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::BubblingGasification<CloudType>::calcXc
(
    const label celli
) const
{
    scalarField Xc(this->owner().thermo().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] =
            this->owner().thermo().carrier().Y()[i][celli]
           /this->owner().thermo().carrier().W(i);
    }

    return Xc/sum(Xc);
}


template<class CloudType>
Foam::scalar Foam::BubblingGasification<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*cbrt(Sc);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BubblingGasification<CloudType>::BubblingGasification
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
    {
        WarningInFunction
            << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        // Determine mapping between liquid and carrier phase species
        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            liqToCarrierMap_[i] =
                owner.composition().carrierId(activeLiquids_[i]);
        }

        // Determine mapping between model active liquids and global liquids
        const label idLiquid = owner.composition().idLiquid();
        forAll(activeLiquids_, i)
        {
            liqToLiqMap_[i] =
                owner.composition().localId(idLiquid, activeLiquids_[i]);
        }
    }
}


template<class CloudType>
Foam::BubblingGasification<CloudType>::BubblingGasification
(
    const BubblingGasification<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BubblingGasification<CloudType>::~BubblingGasification()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BubblingGasification<CloudType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar Re,
    const scalar Pr,
    const scalar d, // Parcel diameter
    const scalar nu,
    const scalar T, // Parcel temperature
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& X,
    scalarField& dMassPC // Returns the parcel mass to gasify
) const
{

//    // immediately evaporate mass that has reached critical condition
//    // prevents temperature of parcel from going above this value
//    // - will release delta spike of mass if this occurs
    scalar Tclip(T);
    
    if ((liquids_.Tc(X) - T) < SMALL)
    {
        if (debug)
        {
            WarningInFunction
                << "Parcel reached critical conditions: "
                << "clip temperature for property evaluations" << endl;
        }

        Tclip = liquids_.Tc(X);

      //  forAll(activeLiquids_, i)
      //  {
      //      const label lid = liqToLiqMap_[i];
      //      dMassPC[lid] = GREAT;
      //  }

      //  return;
    }

//    // droplet surface pressure assumed to surface vapour pressure
//    scalar ps = liquids_.pv(pc, Ts, X);
//
//    // vapour density at droplet surface [kg/m3]
//    scalar rhos = ps*liquids_.W(X)/(RR*Ts);
//
//    // construct carrier phase species volume fractions for cell, celli
//    const scalarField XcMix(calcXc(celli));
//
//    // carrier thermo properties
//    scalar Hsc = 0.0;
//    scalar Hc = 0.0;
//    scalar Cpc = 0.0;
//    scalar kappac = 0.0;
//    forAll(this->owner().thermo().carrier().Y(), i)
//    {
//        scalar Yc = this->owner().thermo().carrier().Y()[i][celli];
//        Hc += Yc*this->owner().thermo().carrier().Ha(i, pc, Tc);
//        Hsc += Yc*this->owner().thermo().carrier().Ha(i, ps, Ts);
//        Cpc += Yc*this->owner().thermo().carrier().Cp(i, ps, Ts);
//        kappac += Yc*this->owner().thermo().carrier().kappa(i, ps, Ts);
//    }

    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];

//        // boiling temperature at cell pressure for liquid species lid [K]
//        const scalar TBoil = liquids_.properties()[lid].pvInvert(pc);
//
//        // limit droplet temperature to boiling/critical temperature
//        const scalar Td = min(T, 0.999*TBoil);
//
//        // saturation pressure for liquid species lid [Pa]
//        const scalar pSat = liquids_.properties()[lid].pv(pc, Td);
//
//        // carrier phase concentration
//        const scalar Xc = XcMix[gid];

        const scalar A(1.0e15); // 1/s
        const scalar Ea(224903.7); // J/mol
        const scalar& Rgas = constant::physicoChemical::R.value(); // J/mol/K 
        const scalar Ta(Ea/Rgas);

        const scalar omega = A*exp(-Ta/T); // 1/s
        const scalar rhod = liquids_.properties()[lid].rho(pc,Tclip);
        const scalar vold = pow(d,3)*pi*(1.0/6.0); 

        dMassPC[lid] += omega*rhod*vold*dt; // has units kg 

//        Info << "parcel evap increment: " << dMassPC << nl
//             << "celli: " << celli << nl
//             << "T: " << T << nl
//             << "rhod: " << rhod << nl
//             << "d: " << d << nl
//             << "vold: " << vold<< nl
//             << "dt: " << dt << nl
//             << "omega (1/s): " << omega << endl;
//        if (Xc*pc > pSat)
//        {
//            // saturated vapour - no phase change
//        }
//        else
//        {
//            // vapour diffusivity [m2/s]
//            const scalar Dab = liquids_.properties()[lid].D(ps, Ts);
//
//            // Schmidt number
//            const scalar Sc = nu/(Dab + ROOTVSMALL);
//
//            // Sherwood number
//            const scalar Sh = this->Sh(Re, Sc);
//
//
//            if (pSat > 0.999*pc)
//            {
//                // boiling
//
//                const scalar deltaT = max(T - TBoil, 0.5);
//
//                // vapour heat of formation
//                const scalar hv = liquids_.properties()[lid].hl(pc, Td);
//
//                // empirical heat transfer coefficient W/m2/K
//                scalar alphaS = 0.0;
//                if (deltaT < 5.0)
//                {
//                    alphaS = 760.0*pow(deltaT, 0.26);
//                }
//                else if (deltaT < 25.0)
//                {
//                    alphaS = 27.0*pow(deltaT, 2.33);
//                }
//                else
//                {
//                    alphaS = 13800.0*pow(deltaT, 0.39);
//                }
//
//                // flash-boil vaporisation rate
//                const scalar Gf = alphaS*deltaT*pi*sqr(d)/hv;
//
//                // model constants
//                // NOTE: using Sherwood number instead of Nusselt number
//                const scalar A = (Hc - Hsc)/hv;
//                const scalar B = pi*kappac/Cpc*d*Sh;
//
//                scalar G = 0.0;
//                if (A > 0.0)
//                {
//                    // heat transfer from the surroundings contributes
//                    // to the vaporisation process
//                    scalar Gr = 1e-5;
//
//                    for (label i=0; i<50; i++)
//                    {
//                        scalar GrDash = Gr;
//
//                        G = B/(1.0 + Gr)*log(1.0 + A*(1.0 + Gr));
//                        Gr = Gf/G;
//
//                        if (mag(Gr - GrDash)/GrDash < 1e-3)
//                        {
//                            break;
//                        }
//                    }
//                }
//
//                dMassPC[lid] += (G + Gf)*dt;
//            }
//            else
//            {
//                // evaporation
//
//                // surface molar fraction - Raoult's Law
//                const scalar Xs = X[lid]*pSat/pc;
//
//                // molar ratio
//                const scalar Xr = (Xs - Xc)/max(SMALL, 1.0 - Xs);
//
//                if (Xr > 0)
//                {
//                    // mass transfer [kg]
//                    dMassPC[lid] += pi*d*Sh*Dab*rhos*log(1.0 + Xr)*dt;
//                }
//            }
//        }
    }
}


template<class CloudType>
Foam::scalar Foam::BubblingGasification<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    scalar dh = 0;

    scalar TDash = T;
    if (liquids_.properties()[idl].pv(p, T) >= 0.999*p)
    {
        TDash = liquids_.properties()[idl].pvInvert(p);
    }

    typedef PhaseChangeModel<CloudType> parent;
    switch (parent::enthalpyTransfer_)
    {
        case (parent::etLatentHeat):
        {
            dh = liquids_.properties()[idl].hl(p, TDash);
            break;
        }
        case (parent::etEnthalpyDifference):
        {
            scalar hc = this->owner().composition().carrier().Ha(idc, p, TDash);
            scalar hp = liquids_.properties()[idl].h(p, TDash);

            dh = hc - hp;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown enthalpyTransfer type" << abort(FatalError);
        }
    }

    return dh;
}


template<class CloudType>
Foam::scalar Foam::BubblingGasification<CloudType>::Tvap
(
    const scalarField& X
) const
{
    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::BubblingGasification<CloudType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return 720.0;// AK - placeholder liquids_.pvInvert(p, X);
}


// ************************************************************************* //
