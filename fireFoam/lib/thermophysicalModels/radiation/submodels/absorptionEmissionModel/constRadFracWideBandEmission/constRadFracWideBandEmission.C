/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "constRadFracWideBandEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"
#include "basicMultiComponentMixture.H"

#include "surfaceFields.H"
#include "radiationModel.H"
#include "fvDOM.H"
#include "blackBodyEmission.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(constRadFracWideBandEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            constRadFracWideBandEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::constRadFracWideBandEmission::constRadFracWideBandEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    radScaling(coeffsDict_.lookupOrDefault<Switch>("radScaling",false)),
    Ehrr1_(coeffsDict_.lookupOrDefault<scalar>("Ehrr1",0.3)),
    Ehrr2_
    (
        coeffsDict_.lookupOrDefault<scalar>("Ehrr2",0.3)
    ),
    patchName1_(coeffsDict_.lookup("patch1")),
    patchName2_(coeffsDict_.lookup("patch2")),
    EhrrMin_
    (
        min(Ehrr1_, Ehrr2_)
    ),
    EhrrMax_
    (
        max(Ehrr1_, Ehrr2_)
    ),
    radTemp_(coeffsDict_.lookupOrDefault<scalar>("RadTemperature",1400.)),
    fracsSet_(false)
{
    // Read band limits from multimedia-radiation file
    const IOdictionary propDict
    (
	IOobject
	(
	    "multimediaRadProperties",
	    this->mesh().time().constant(),
	    this->mesh(), 
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	)
    );
    (propDict.subDict("allBands")).lookup("bandLimits") >> iBands_;
    nBands_ = iBands_.size();

    enFracs_.setSize(nBands_);


}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::constRadFracWideBandEmission::~constRadFracWideBandEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::constRadFracWideBandEmission::aCont(const label bandI) const
{

    tmp<volScalarField> a
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    return a;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::constRadFracWideBandEmission::eCont(const label bandI) const
{
    tmp<volScalarField> e
    (
        new volScalarField
        (
            IOobject
            (
                "eCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    return e;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::constRadFracWideBandEmission::ECont(const label bandI) const
{

    if (!fracsSet_)
    {
     Foam::radiation::constRadFracWideBandEmission& modAbsModelObj = const_cast<Foam::radiation::constRadFracWideBandEmission&>(*this);
     modAbsModelObj.setEnergyFracs();
//        FatalErrorInFunction
//            << "Energy fractions need to be set by this point. Error in the logic. Please check." 
//            << abort(FatalError); 
    }

    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    scalar RadFraction = 0;

    if (radScaling)
    {

        //TODO: this doesn't need to be recomputed for each ILambda solve

        const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");

        scalar mlr1(0.0);

        // kvm added this ....

        forAll(patchName1_,i)
        {
            const label patchI = mesh_.boundaryMesh().findPatchID(patchName1_[i]);
            if(patchI<0)
            {
                FatalErrorIn("radScaling.H")
                    << "patch " << patchName1_[i] << " not found" << nl
                    << abort(FatalError);
            }
            mlr1 += -gSum(phi.boundaryField()[patchI]);
        }

        scalar mlr2(0.0);
        scalar mlr2xEhrr2(0.0);

        forAll(patchName2_,i)
        {
            const label patchI = mesh_.boundaryMesh().findPatchID(patchName2_[i]);
            if(patchI<0)
            {
                FatalErrorIn("radScaling.H")
                    << "patch " << patchName2_[i] << " not found" << nl
                    << abort(FatalError);
            }

            {
                scalar mlr2p = -gSum(phi.boundaryField()[patchI]);
                mlr2 += mlr2p;
                mlr2xEhrr2 += mlr2p*Ehrr2_;
            }

        }

        if(debug)
        {
            Info << "mlr for patches " << patchName1_ << " is " << mlr1 << endl;
            Info << "mlr for patches " << patchName2_ << " is " << mlr2 << endl;
        }

        RadFraction = (mlr1*Ehrr1_ + mlr2xEhrr2)
                    / max(SMALL, (mlr1 + mlr2));
        RadFraction = max(EhrrMin_,RadFraction);

    }
    else
    {
        RadFraction = EhrrCoeff_;
    }


    if (mesh_.foundObject<volScalarField>("Qdot"))
    {
        const volScalarField& Qdot =
            mesh_.lookupObject<volScalarField>("Qdot");
        if (Qdot.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E.ref().ref() = RadFraction*Qdot*enFracs_[bandI];
        }
        else
        {
            Info << "Qdot dimensions incorrect" << endl;
        }

        static word timeName = "null";
        if (timeName != mesh().time().timeName())
        {
            Info << "Radiant Fraction is " << RadFraction << endl;
            timeName = mesh().time().timeName();
        }
    }

    return E;
}


void Foam::radiation::constRadFracWideBandEmission::correct
(
    Foam::volScalarField& a,
    PtrList<Foam::volScalarField>& aLambda
) 
{

    a = dimensionedScalar("zero", dimless/dimLength, 0.0);

    if (!fracsSet_)
    {
      Foam::radiation::constRadFracWideBandEmission::setEnergyFracs();  // Setting energy fractions for various bands 
      // setEnergyFracs();
    }

}


void Foam::radiation::constRadFracWideBandEmission::setEnergyFracs()
{
    const radiationModel& radiation =
        mesh().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));
     
    const scalar radTemp = radTemp_;
    for (label j=0; j<nBands_; j++)
    {
       const Vector2D<scalar>& iBand = iBands_[j];

       enFracs_[j] = dom.blackBody().fDeltaLambdaT(radTemp,iBand); 
    } 
    
    fracsSet_ = true; 

}
// ************************************************************************* //
