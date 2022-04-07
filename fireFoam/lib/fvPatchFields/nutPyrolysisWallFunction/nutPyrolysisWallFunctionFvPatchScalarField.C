/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "nutPyrolysisWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"

#define DEBUG(x) {                                              \
            std::streamsize p = std::cout.precision();              \
            std::ios::fmtflags myFlags;                             \
            myFlags = cout.flags();                                 \
            std::cout.precision(10);                                \
            std::cout.setf(std::ios::fixed,std::ios::floatfield);   \
            std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
            std::cout << "p" << Pstream::myProcNo();                \
            std::cout << " " << #x " = " << x << std::endl;         \
            std::cout.precision(p);                                 \
            std::cout.flags(myFlags);                               \
        }
#define TRACE(s) {                                              \
            std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
            std::cout << "p" << Pstream::myProcNo();                \
            std::cout << " " << #s << std::endl;                    \
            s;                                                      \
        }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutPyrolysisWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const compressible::turbulenceModel& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& rhow = patch().lookupPatchField<volScalarField, scalar>("rho");

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField& yw = turbModel.y()[patchi];

    const scalarField& phiw = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalarField pyrolysisMassFlux( - fuelConversionRatio_*phiw/patch().magSf()); // ConvertionRatio=2.5 for bio-fuel

    //scalarField nutSF = db().template lookupObject<volScalarField>("nut");
    volScalarField nutSF = db().lookupObject<volScalarField>("nut");
    //const labelUList& faceCells = nutSF.boundaryField()[patchi].faceCells();//patch().faceCells()[patchi];
    const labelUList& faceCells = patch().faceCells();//patch().faceCells()[patchi];

    const vectorField& faceC = patch().Cf();

    //volScalarField PV = db().lookupObject<volScalarField>("PV");

    forAll(nutw,facei)
    {
        scalar uyplus0(1.0);
        scalar U0(magUp[facei]+1.0e-2);
        scalar Y0(yw[facei]);
        scalar rho0(rhow[facei]);
        //scalar nu0(nuwInt[facei]);
        scalar nu0(nuw[facei]);
        scalar mu0(nu0*rho0);
        scalar mf0(max(1.0e-8,pyrolysisMassFlux[facei]));
        //scalar mf0(1.0e-8);
        scalar tau0(mu0*U0/Y0);

        scalar Ustar0(sqrt(tau0/rho0));
        scalar Mstar0(sqrt(tau0*rho0));
        scalar Uplus0(U0/Ustar0);
        scalar Yplus0(Y0*Ustar0/nu0);
        scalar Mplus0(mf0/Mstar0);

        scalar YplusErr(10);
        scalar YplusOld(1);
        scalar kEy(1.0);
        int j(0);
        while((j<10) && (YplusErr>5))
        {
            //Info<<"faceI: "<<facei<<tab<<j<<tab<<YplusErr<<endl;
            //Info<<Yplus0<<tab<<Uplus0<<tab<<Ustar0<<tab<<Mstar0<<tab<<mf0<<tab<<U0<<tab<<Y0<<tab<<tau0<<endl;
            if(Yplus0 > 12)
            {
                kEy = log(E_*Yplus0)/kappa_;
                Uplus0 = kEy + Mplus0*kEy*kEy/Cb_;
                Ustar0 = U0/Uplus0;
                Yplus0 = Y0*Ustar0/nu0;
                tau0   = Ustar0*Ustar0*rho0;
                Mstar0 = sqrt(tau0*rho0);
                Mplus0 = mf0/Mstar0;
            }
            else
            {
                Uplus0 = (exp(mf0*Y0/mu0)-1.0)/(Mplus0 + 1.0e-6);
                Ustar0 = U0/Uplus0;
                Yplus0 = Y0*Ustar0/nu0;
                tau0   = U0*mf0/(exp(mf0*Y0/mu0) - 1.0 + 1.0e-6);
                Mstar0 = sqrt(tau0*rho0);
                Mplus0 = mf0/Mstar0;
            }
            j++;
            YplusErr = Yplus0 - YplusOld;
            YplusOld = Yplus0;
        }

        uyplus0 = Yplus0/Uplus0;

        nutw[facei] = nuw[facei]*max(0,(uyplus0 - 1.0));
        //Info<<"dbg--> "<<facei<<tab<<j<<tab<<Yplus0<<tab<<Uplus0<<tab<<nutw[facei]<<tab<<mf0<<endl;
        //Info<<"dbg--> "<<mf0*1000<<tab<<Yplus0<<tab<<Uplus0<<tab<<nutw[facei]
        //               <<tab<<nu0<<tab<<U0<<endl;
                       //<<tab<<nu0<<tab<<mu0<<tab<<rho0<<tab<<U0<<endl;

        //nutw[facei] = nuw[facei]*(UYplus(magUp[facei],y0[facei],fuelMassFlux[facei]) - 1.0);

        //- Temporary test (set nuw = 1.2e-4*BlowingEffect)
        //scalar mRatio(fuelMassFlux[facei]/0.010 + 0.001);
        //nutw[facei] = 1.4e-4*mRatio/(exp(mRatio)-1.0);
        //nutw[facei] = nutw[facei] + nuw[facei]/(PV[faceCells[facei]] + SMALL);

        if(correctInternal_)
        {
            nutSF[faceCells[facei]] = nutw[facei];
        }

        //if((faceC[facei].x() >= 0.0) && (faceC[facei].x() < 0.03))
        //{
        //    //Info<<"dbg--> "<<mf0*1000<<tab<<Yplus0<<tab<<Uplus0<<tab<<nutw[facei]
        //    //               <<tab<<nu0<<tab<<U0<<endl;
        //    Info<<"dbg--> "<<facei<<tab<<faceC[facei].z()<<tab<<PV[faceCells[facei]]<<tab<<nutw[facei]<<endl;
        //}
    }



    //forAll(yPlus, facei)
    //{
    //    if (yPlus[facei] > yPlusLam_)
    //    {
    //        nutw[facei] = nuw[facei]*(yPlus[facei]/UPlus(yPlus[facei]) - 1.0);
    //    }
    //}

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutPyrolysisWallFunctionFvPatchScalarField::nutPyrolysisWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    Cb_(2.5),
    fuelConversionRatio_(2.5),
    correctInternal_(false)
{}


nutPyrolysisWallFunctionFvPatchScalarField::nutPyrolysisWallFunctionFvPatchScalarField
(
    const nutPyrolysisWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Cb_(ptf.Cb_),
    fuelConversionRatio_(ptf.fuelConversionRatio_),
    correctInternal_(ptf.correctInternal_)
{}


nutPyrolysisWallFunctionFvPatchScalarField::nutPyrolysisWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    Cb_(dict.lookupOrDefault<scalar>("Cb", 2.5)),
    fuelConversionRatio_(dict.lookupOrDefault<scalar>("fuelConversionRatio", 2.5)),
    correctInternal_(dict.lookupOrDefault<bool>("correctInternal", false))
{}


nutPyrolysisWallFunctionFvPatchScalarField::nutPyrolysisWallFunctionFvPatchScalarField
(
    const nutPyrolysisWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    Cb_(wfpsf.Cb_),
    fuelConversionRatio_(wfpsf.fuelConversionRatio_),
    correctInternal_(wfpsf.correctInternal_)
{}


nutPyrolysisWallFunctionFvPatchScalarField:: nutPyrolysisWallFunctionFvPatchScalarField
(
    const nutPyrolysisWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    Cb_(wfpsf.Cb_),
    fuelConversionRatio_(wfpsf.fuelConversionRatio_),
    correctInternal_(wfpsf.correctInternal_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutPyrolysisWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}

//tmp<scalarField> nutPyrolysisWallFunctionFvPatchScalarField::UYplus() const
//{
//    const label patchi = patch().index();
//
//    const compressible::turbulenceModel& turbModel = db().lookupObject<compressible::turbulenceModel>
//    (
//        IOobject::groupName
//        (
//            turbulenceModel::propertiesName,
//            internalField().group()
//        )
//    );
//    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
//    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
//    const scalarField& y0 = turbModel.y()[patchi];
//
//    const scalarField& phiw = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
//    scalarField fuelMassFlux( - phiw/patch().magSf()); //convert to g/m2/s, and back to pyrolysate
//
//    return y0;
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutPyrolysisWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
