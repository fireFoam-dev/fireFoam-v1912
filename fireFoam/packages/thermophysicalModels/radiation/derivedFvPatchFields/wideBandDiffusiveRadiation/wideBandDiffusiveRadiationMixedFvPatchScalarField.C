/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "wideBandDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
#include "wideBandAbsorptionEmission.H"
#include "constants.H"
#include "boundaryRadiationProperties.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
wideBandDiffusiveRadiationMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    const label patchi = patch().index();

    if (dom.nLambda() == 0)
    {
        FatalErrorInFunction
            << " a non-grey boundary condition is used with a grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;
    const vectorField n(patch().Sf()/patch().magSf());

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    const scalarField nAve(n & ray.dAve());

    // ankur
    ray.qrLambda(lambdaId).boundaryFieldRef()[patchi] += Iw*nAve;

    // ankur
    const scalarField pT
    (
        dom.T().boundaryField()[patchi]
    );

    const scalarField Eb
    (
        dom.enFracLambda(lambdaId).boundaryField()[patchi]*Foam::pow(pT,4.)*physicoChemical::sigma.value()  // ankur
        //dom.blackBody().bLambda(lambdaId).boundaryField()[patchi] // making code generic for wide-band and WSGG banded models.. 
    );

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(internalField().mesh());


    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patch().index(), lambdaId)
    );

    const scalarField& emissivity = temissivity();

    const tmp<scalarField> ttransmissivity
    (
        boundaryRadiation.transmissivity(patch().index(), lambdaId)
    );
    const scalarField& transmissivity = ttransmissivity();

    scalarField& qem = ray.qemLambda(lambdaId).boundaryFieldRef()[patchi]; // ankur
    scalarField& qin = ray.qinLambda(lambdaId).boundaryFieldRef()[patchi]; // ankur

    // Use updated Ir while iterating over rays
    // avoids to used lagged qin
    scalarField Ir = dom.IRay(0).qinLambda(lambdaId).boundaryField()[patchi]; // ankur

    for (label rayI=1; rayI < dom.nRay(); rayI++)
    {
        Ir += dom.IRay(rayI).qinLambda(lambdaId).boundaryField()[patchi]; // ankur
    }
    // kvm, look in to this block and see if we need it to replace the block above, 2018-04
    // Calculate Ir into the wall on the same lambdaId
    // scalarField Ir(patch().size(), 0.0);
    // forAll(Iw, facei)
    // {
    //     for (label rayi=0; rayi < dom.nRay(); rayi++)
    //     {
    //         const vector& d = dom.IRay(rayi).d();

    //         if ((-n[facei] & d) < 0.0)
    //         {
    //             // q into the wall
    //             const scalarField& IFace =
    //                 dom.IRay(rayi).ILambda(lambdaId).boundaryField()[patchi];

    //             const vector& rayDave = dom.IRay(rayi).dAve();
    //             Ir[facei] += IFace[facei]*(n[facei] & rayDave);
    //         }
    //     }
    // }

    if (dom.useSolarLoad())
    {
        // Looking for primary heat flux single band
        Ir += patch().lookupPatchField<volScalarField,scalar>
        (
            dom.primaryFluxName_ + "_" + name(lambdaId)
        );

        word qSecName = dom.relfectedFluxName_ + "_" + name(lambdaId);

        if (this->db().foundObject<volScalarField>(qSecName))
        {
             const volScalarField& qSec =
                this->db().lookupObject<volScalarField>(qSecName);

            Ir += qSec.boundaryField()[patch().index()];
        }
    }

    scalarField Iexternal(this->size(), 0.0);

    if (dom.useExternalBeam())
    {
        const vector sunDir = dom.solarCalc().direction();
        const scalar directSolarRad =
            dom.solarCalc().directSolarRad()
           *dom.spectralDistribution()[lambdaId];

        //label nRaysBeam = dom.nRaysBeam();
        label SunRayId(-1);
        scalar maxSunRay = -GREAT;

        // Looking for the ray closest to the Sun direction
        for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
            const vector& iD = dom.IRay(rayI).d();
            scalar dir = sunDir & iD;
            if (dir > maxSunRay)
            {
                maxSunRay = dir;
                SunRayId = rayI;
            }
        }

        if (rayId == SunRayId)
        {
            const scalarField nAve(n & dom.IRay(rayId).dAve());
            forAll(Iexternal, faceI)
            {
                Iexternal[faceI] = directSolarRad/mag(dom.IRay(rayId).dAve());
            }
        }
    }

    forAll(Iw, facei)
    {
        const vector& d = dom.IRay(rayId).d();

        if ((-n[facei] & d) > 0.0)
        {
            // direction out of the wall
            refGrad()[facei] = 0.0;
            valueFraction()[facei] = 1.0;
            refValue()[facei] =
                Iexternal[facei]*transmissivity[facei]
              + (
                    Ir[facei]*(1.0 - emissivity[facei])
                  + emissivity[facei]*Eb[facei]
                )/pi;

            // Emitted heat flux from this ray direction (sum over lambdaId)
            qem[facei] += refValue()[facei]*nAve[facei];
        }
        else
        {
            // direction into the wall
            valueFraction()[facei] = 0.0;
            refGrad()[facei] = 0.0;
            refValue()[facei] = 0.0; //not used

            // Incident heat flux on this ray direction (sum over lambdaId)
            qin[facei] += Iw[facei]*nAve[facei];
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wideBandDiffusiveRadiationMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
