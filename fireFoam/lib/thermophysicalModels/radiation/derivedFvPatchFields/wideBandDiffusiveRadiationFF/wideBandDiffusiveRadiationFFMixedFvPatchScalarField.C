/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "wideBandDiffusiveRadiationFFMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
//#include "wideBandAbsorptionEmission.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::
wideBandDiffusiveRadiationFFMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBaseFF(p, "undefined", "undefined", scalarField::null(), scalarField::null()),
    TName_("T")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::
wideBandDiffusiveRadiationFFMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationFFMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBaseFF
    (
        p,
        ptf.emissivityMethod(),
        ptf.absorptivityMethod(),
        ptf.emissivity_,
        ptf.absorptivity_
    ),
    TName_(ptf.TName_)
{}


Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::
wideBandDiffusiveRadiationFFMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBaseFF(p, dict),
    TName_(dict.lookupOrDefault<word>("T", "T"))
{
    if (dict.found("value"))
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
        const scalarField& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        refValue() =
            4.0*physicoChemical::sigma.value()*pow4(Tp)*emissivity()/pi;
        refGrad() = 0.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::
wideBandDiffusiveRadiationFFMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationFFMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    radiationCoupledBaseFF
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.absorptivityMethod(),
        ptf.emissivity_,
        ptf.absorptivity_
    ),
    TName_(ptf.TName_)
{}


Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::
wideBandDiffusiveRadiationFFMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationFFMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    radiationCoupledBaseFF
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.absorptivityMethod(),
        ptf.emissivity_,
        ptf.absorptivity_
    ),
    TName_(ptf.TName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::
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

    ray.qr().boundaryFieldRef()[patchi] += Iw*nAve;

    const scalarField pT
    (
        dom.T().boundaryField()[patchi]
    );

    const scalarField Eb
    (
        dom.enFracLambda(lambdaId).boundaryField()[patchi]*Foam::pow(pT,4.)*physicoChemical::sigma.value()
        //dom.blackBody().bLambda(lambdaId).boundaryField()[patchi]
    );

    scalarField temissivity = emissivity();

    scalarField tabsorptivity;

    if (absMethod_ == EMISSIVITY)
    {
    tabsorptivity = temissivity;
    }
    else
    {
    tabsorptivity = absorptivity();
    }

    scalarField& qem = ray.qem().boundaryFieldRef()[patchi];
    scalarField& qin = ray.qin().boundaryFieldRef()[patchi];

    // Use updated Ir while iterating over rays
    // avoids to used lagged qin
    scalarField Ir = dom.IRay(0).qin().boundaryField()[patchi];

    for (label rayI=1; rayI < dom.nRay(); rayI++)
    {
        Ir += dom.IRay(rayI).qin().boundaryField()[patchi];
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
                (
                    Ir[facei]*(1.0 - tabsorptivity[facei])
                  + temissivity[facei]*Eb[facei]
                )/pi;

            // Emmited heat flux from this ray direction
            qem[facei] = refValue()[facei]*nAve[facei];
        }
        else
        {
            // direction into the wall
            valueFraction()[facei] = 0.0;
            refGrad()[facei] = 0.0;
            refValue()[facei] = 0.0; //not used

            // Incident heat flux on this ray direction
            qin[facei] = Iw[facei]*nAve[facei];
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::wideBandDiffusiveRadiationFFMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    radiationCoupledBaseFF::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wideBandDiffusiveRadiationFFMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
