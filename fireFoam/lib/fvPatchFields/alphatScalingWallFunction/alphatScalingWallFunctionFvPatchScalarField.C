/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "alphatScalingWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
//#include "muSgsBuoyantWallFunctionFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

/*
scalar alphatScalingWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphatScalingWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatScalingWallFunctionFvPatchScalarField::maxIters_ = 10;
*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*
void alphatScalingWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphatScalingWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatScalingWallFunctionFvPatchScalarField::
alphatScalingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    floorSurface_(0.6),
    fuelConversionRatio_(2.5),
    m60_(0.01),
    Prt_(1.0),
    correctInternal_(false),
    Cdelta1_(175.0),
    Cdelta2_(800.0),
    deltaMin_(0.0015),
    delta1_(0.0015),
    delta2_(0.01)
    //Gc_
    //(
    //    IOobject
    //    (
    //        "Gc",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //),
    //Bc_
    //(
    //    IOobject
    //    (
    //        "Bc",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //),
    //Cc_
    //(
    //    IOobject
    //    (
    //        "Cc",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //),
    //Cc2_
    //(
    //    IOobject
    //    (
    //        "Cc2",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //)
{
    //checkType();
    //read();
}


alphatScalingWallFunctionFvPatchScalarField::
alphatScalingWallFunctionFvPatchScalarField
(
    const alphatScalingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    floorSurface_(ptf.floorSurface_),
    fuelConversionRatio_(ptf.fuelConversionRatio_),
    m60_(ptf.m60_),
    Prt_(ptf.Prt_),
    correctInternal_(ptf.correctInternal_),
    minCombustionCorr_(ptf.minCombustionCorr_),
    Cdelta1_(ptf.Cdelta1_),
    Cdelta2_(ptf.Cdelta2_),
    deltaMin_(ptf.deltaMin_),
    delta1_(ptf.delta1_),
    delta2_(ptf.delta2_)
    //Gc_(ptf.Gc_),
    //Bc_(ptf.Bc_),
    //Cc_(ptf.Cc_),
    //Cc2_(ptf.Cc2_)
{}


alphatScalingWallFunctionFvPatchScalarField::
alphatScalingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    floorSurface_(dict.lookupOrDefault<scalar>("floorSurface", 0.6)),
    fuelConversionRatio_(dict.lookupOrDefault<scalar>("fuelConversionRatio", 2.5)),
    m60_(dict.lookupOrDefault<scalar>("m60", 0.01)),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 1.0)),
    correctInternal_(dict.lookupOrDefault<bool>("correctInternal", false)),
    minCombustionCorr_(dict.lookupOrDefault<bool>("minCombustionCorr", false)),
    Cdelta1_(dict.lookupOrDefault<scalar>("Cdelta1", 175.0)),
    Cdelta2_(dict.lookupOrDefault<scalar>("Cdelta2", 800.0)),
    deltaMin_(dict.lookupOrDefault<scalar>("deltaMin", 0.0015)),
    delta1_(dict.lookupOrDefault<scalar>("delta1", 0.0015)),
    delta2_(dict.lookupOrDefault<scalar>("delta2", 0.01))
    //Gc_
    //(
    //    IOobject
    //    (
    //        "Gc",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //),
    //Bc_
    //(
    //    IOobject
    //    (
    //        "Bc",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //),
    //Cc_
    //(
    //    IOobject
    //    (
    //        "Cc",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //),
    //Cc2_
    //(
    //    IOobject
    //    (
    //        "Cc2",
    //        this->mesh().time().timeName(),
    //        this->mesh(),
    //        IOobject::NO_READ,
    //        //IOobject::AUTO_WRITE
    //        IOobject::NO_WRITE
    //    ),
    //    this->mesh(),
    //    dimensionedScalar("zero", dimless, 1.0)
    //)
{
//    checkType();
    //read();
}


alphatScalingWallFunctionFvPatchScalarField::
alphatScalingWallFunctionFvPatchScalarField
(
    const alphatScalingWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    floorSurface_(tppsf.floorSurface_),
    fuelConversionRatio_(tppsf.fuelConversionRatio_),
    m60_(tppsf.m60_),
    Prt_(tppsf.Prt_),
    correctInternal_(tppsf.correctInternal_),
    minCombustionCorr_(tppsf.minCombustionCorr_),
    Cdelta1_(tppsf.Cdelta1_),
    Cdelta2_(tppsf.Cdelta2_),
    deltaMin_(tppsf.deltaMin_),
    delta1_(tppsf.delta1_),
    delta2_(tppsf.delta2_)
{
//    checkType();
}


alphatScalingWallFunctionFvPatchScalarField::
alphatScalingWallFunctionFvPatchScalarField
(
    const alphatScalingWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    floorSurface_(tppsf.floorSurface_),
    fuelConversionRatio_(tppsf.fuelConversionRatio_),
    m60_(tppsf.m60_),
    Prt_(tppsf.Prt_),
    correctInternal_(tppsf.correctInternal_),
    minCombustionCorr_(tppsf.minCombustionCorr_),
    Cdelta1_(tppsf.Cdelta1_),
    Cdelta2_(tppsf.Cdelta2_),
    deltaMin_(tppsf.deltaMin_),
    delta1_(tppsf.delta1_),
    delta2_(tppsf.delta2_)
    //Gc_(tppsf.Gc_),
    //Bc_(tppsf.Bc_),
    //Cc_(tppsf.Cc_),
    //Cc2_(tppsf.Cc2_)
{
//    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatScalingWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
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

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );

    const fvPatchScalarField& alphaw = turbModel.alpha()().boundaryField()[patchi];
    scalarField& alphatw = *this;

    //const scalarField& nutw = patch().lookupPatchField<volScalarField, scalar>("nut");
    //const scalarField& rhow = patch().lookupPatchField<volScalarField, scalar>("rho");
    //const scalarField Cpw = thermo.Cp()().boundaryField()[patchi];
    const fvPatchScalarField& Tw = thermo.T().boundaryField()[patchi];
    const scalarField T(Tw.patchInternalField());
    //const scalarField magGradTw(max(mag(Tw.snGrad()), VSMALL));
    //const scalarField gradTw(Tw.snGrad());

    volScalarField alphatSF = db().lookupObject<volScalarField>("alphat");
    const labelUList& faceCells = patch().faceCells();
    const scalarField& yw = turbModel.y()[patchi];

    const vectorField& faceC = patch().Cf();
    const vectorField& faceSf = patch().Sf();
    const scalarField& faceSfMag = patch().magSf();
    //vector verticalDirection(0 0 1);
    volScalarField PV = db().lookupObject<volScalarField>("PV");


    const scalarField& phiw = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalarField pyrolysisMassFlux( - fuelConversionRatio_*phiw/patch().magSf());
    //scalarField& gridC = patch().lookupPatchField<surfaceScalarField, scalar>("GridC");
    //scalarField& blowC = patch().lookupPatchField<surfaceScalarField, scalar>("BlowC");
    //scalarField& combC = patch().lookupPatchField<surfaceScalarField, scalar>("CombC");
    //scalarField& combC2 = patch().lookupPatchField<surfaceScalarField, scalar>("CombC2");

    forAll(alphatw,facei)
    {
        //alphatw[facei]=Prt_*nutw[facei]*rhow[facei];
        scalar mRatio(pyrolysisMassFlux[facei]/m60_ + 0.001);
        scalar BF(mRatio/(exp(mRatio)-1.0));
        scalar PVmin(max(0.1,(PV[faceCells[facei]]+SMALL)));
        scalar dDelta(max(0,yw[facei]*2.0-deltaMin_));
        scalar dDelta1(max(0,yw[facei]*2.0-delta1_));
        scalar dDelta2(delta2_- delta1_);

        scalar Gc(Cdelta1_*dDelta - Cdelta2_*dDelta*dDelta);
        //scalar Bc(min(1.0,max(BF,(1.0-dDelta1*(1.0-BF)/dDelta2))));
        scalar Bc(BF + (1-BF)*exp(-dDelta1/dDelta2));

        //scalar ExtraCorr(min(20,(T[facei]-300*(1.0-PVmin)-Tw[facei]*PVmin)/max(1.0,(T[facei]-Tw[facei]))));
        scalar ExtraCorr(1.0);
        scalar Tcorr(300.0 + (T[facei] - 300.0)/PVmin);
        //if((fabs(T[facei]-Tw[facei])<1.0) && (T[facei] > 400))
        if(fabs(T[facei]-Tw[facei])<1.0)
        {
            //Info<<"Surface T = Tinternal "<<facei<<tab<<T[facei]<<tab<<Tw[facei]<<tab<<Tcorr<<endl;
            ExtraCorr = 1.0;
        }
        else
        {
            ExtraCorr  = min(10,PVmin*(Tcorr - Tw[facei])/(T[facei] - Tw[facei]));
            if(minCombustionCorr_)
            {
                ExtraCorr = max(0.0, ExtraCorr);
            }
            //Info<<"dbg-facei "<<facei<<tab<<Tw[facei]<<tab<<T[facei]<<tab<<Tcorr<<tab<<PVmin<<tab<<ExtraCorr<<endl;
        }
        //ExtraCorr = 1;

        //Info<<"location: "<<facei<<tab<<faceC[facei].x()<<tab<<faceC[facei].y()<<tab<<faceC[facei].z()<<endl;
        //if(fabs(faceC[facei].x() == 0) && (faceC[facei].y() == 0))
        //{
        //    Info<<"horizontal surface: "<<facei<<tab<<faceC[facei].z()<<endl;
        //}
        scalar horizontalFactor(1.0);
        if(fabs(faceSf[facei].z()/faceSfMag[facei] + 1.0) < 1e-3)
        {
            //Info<<"horizontal surface: "<<facei<<tab<<faceSf[facei].x()<<tab<<faceSf[facei].y()<<tab<<faceSf[facei].z()<<endl;
            horizontalFactor = floorSurface_;
        }
        //Info<<"vector: "<<facei<<tab<<faceSf[facei].x()<<tab<<faceSf[facei].y()<<tab<<faceSf[facei].z()<<tab<<faceSfMag[facei]<<endl;

        alphatw[facei] = horizontalFactor*alphaw[facei]*(Bc*(Gc + 1.0)*ExtraCorr/PVmin - 1.0);

        //gridC[facei] = Gc;
        //blowC[facei] = Bc;
        //combC[facei] = 1.0/PVmin;

        //alphatw[facei] = alphaw[facei]*(BF*(Cdelta_*(yw[facei]*2.0-deltaMin_)+1.0)/PVmin - 1.0);

        if(correctInternal_)
        {
            alphatSF[faceCells[facei]]=alphatw[facei];
        }

    }
}

void alphatScalingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    //writeLocalEntries(os);
    os.writeKeyword("floorSurface") << floorSurface_ << token::END_STATEMENT << nl;
    os.writeKeyword("fuelConversionRatio") << fuelConversionRatio_ << token::END_STATEMENT << nl;
    os.writeKeyword("m60") << m60_ << token::END_STATEMENT << nl;
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
    //os.writeKeyword("correctInternal") << correctInternal_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cdelta1") << Cdelta1_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cdelta2") << Cdelta2_ << token::END_STATEMENT << nl;
    os.writeKeyword("deltaMin") << deltaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("delta1") << delta1_ << token::END_STATEMENT << nl;
    os.writeKeyword("delta2") << delta2_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatScalingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
