/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "greyDiffusiveRadiationFFMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::
greyDiffusiveRadiationFFMixedFvPatchScalarField
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


Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::
greyDiffusiveRadiationFFMixedFvPatchScalarField
(
    const greyDiffusiveRadiationFFMixedFvPatchScalarField& ptf,
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


Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::
greyDiffusiveRadiationFFMixedFvPatchScalarField
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


Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::
greyDiffusiveRadiationFFMixedFvPatchScalarField
(
    const greyDiffusiveRadiationFFMixedFvPatchScalarField& ptf
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


Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::
greyDiffusiveRadiationFFMixedFvPatchScalarField
(
    const greyDiffusiveRadiationFFMixedFvPatchScalarField& ptf,
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

void Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::
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

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    const label patchI = patch().index();

    if (dom.nLambda() != 1)
    {
        FatalErrorIn
        (
            "Foam::radiation::"
            "greyDiffusiveRadiationFFMixedFvPatchScalarField::updateCoeffs"
        )   << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;
    const vectorField n(patch().nf());

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    const scalarField nAve(n & ray.dAve());

    ray.qr().boundaryFieldRef()[patchI] += Iw*nAve;

    const scalarField temissivity = emissivity();

    const scalarField fTransmissivity = filmTransmissivity();
  
    scalarField tabsorptivity;
 
    if (absMethod_ == EMISSIVITY)
    {
    tabsorptivity = temissivity;  
    }
    else
    { 
    tabsorptivity = absorptivity();
    }

    scalarField& qem = ray.qem().boundaryFieldRef()[patchI];
    scalarField& qin = ray.qin().boundaryFieldRef()[patchI];

    const vector& myRayId = dom.IRay(rayId).d();

    // Use updated Ir while iterating over rays
    // avoids to used lagged qin
    scalarField Ir = dom.IRay(0).qin().boundaryField()[patchI];

    for (label rayI=1; rayI < dom.nRay(); rayI++)
    {
        Ir += dom.IRay(rayI).qin().boundaryField()[patchI];
    }

    forAll(Iw, faceI)
    {
        if ((-n[faceI] & myRayId) > 0.0)
        {
            // direction out of the wall
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            //- luwi-- total emissive intensity from solid boundary
            scalar Iem_faceI =
                (
                    fTransmissivity[faceI]*Ir[faceI]*(scalar(1.0) - tabsorptivity[faceI])
                  + temissivity[faceI]*physicoChemical::sigma.value()
                  * pow4(Tp[faceI])
                )/pi;

            // --luwi-- boundary source for RTE is net emission through water film covering solid surface
            refValue()[faceI] = Iem_faceI*fTransmissivity[faceI];
            // Emmited heat flux from this ray direction
            // --luwi-- stored qem is total emission from solid surface, not from its covering water film
            qem[faceI] = Iem_faceI*nAve[faceI];
        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used

            // Incident heat flux on this ray direction
            qin[faceI] = Iw[faceI]*nAve[faceI];
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::greyDiffusiveRadiationFFMixedFvPatchScalarField::write
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
        greyDiffusiveRadiationFFMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
