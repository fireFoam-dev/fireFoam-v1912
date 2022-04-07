/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "wavyFilmHeightFvPatchScalarField.H"
#include "volFields.H"
#include "singleLayerRegion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wavyFilmHeightFvPatchScalarField::wavyFilmHeightFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    mean_(),
    a_(),
    omega_()
{}


Foam::wavyFilmHeightFvPatchScalarField::wavyFilmHeightFvPatchScalarField
(
    const wavyFilmHeightFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    mean_(ptf.mean_().clone().ptr()),
    a_(ptf.a_().clone().ptr()),
    omega_(ptf.omega_().clone().ptr())
{}


Foam::wavyFilmHeightFvPatchScalarField::wavyFilmHeightFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    mean_(DataEntry<scalar>::New("mean", dict)),
    a_(DataEntry<scalar>::New("a", dict)),
    omega_(DataEntry<scalar>::New("omega", dict))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        const scalar t = db().time().timeOutputValue();
        fvPatchScalarField::operator=(mean_->value(t));
    }
}


Foam::wavyFilmHeightFvPatchScalarField::wavyFilmHeightFvPatchScalarField
(
    const wavyFilmHeightFvPatchScalarField& wfhpsf
)
:
    fixedValueFvPatchScalarField(wfhpsf),
    mean_(wfhpsf.mean_().clone().ptr()),
    a_(wfhpsf.a_().clone().ptr()),
    omega_(wfhpsf.omega_().clone().ptr())
{}


Foam::wavyFilmHeightFvPatchScalarField::wavyFilmHeightFvPatchScalarField
(
    const wavyFilmHeightFvPatchScalarField& wfhpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfhpsf, iF),
    mean_(wfhpsf.mean_().clone().ptr()),
    a_(wfhpsf.a_().clone().ptr()),
    omega_(wfhpsf.omega_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wavyFilmHeightFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const scalar t = db().time().timeOutputValue();

    // retrieve the film region from the database

    const regionModels::regionModel& region =
        db().time().lookupObject<regionModels::regionModel>
        (
            "surfaceFilmProperties"
        );

    const regionModels::singleLayerRegion& film =
        dynamic_cast<const regionModels::singleLayerRegion&>(region);

    // calculate the vector tangential to the patch

    const volVectorField& nHat = film.nHat();

    const vectorField nHatp(nHat.boundaryField()[patchI].patchInternalField());

    vectorField nTan(nHatp ^ patch().nf());
    nTan /= mag(nTan) + ROOTVSMALL;

    // calculate distance in patch tangential direction

    const vectorField& Cf = patch().Cf();
    scalarField d(nTan & Cf);
    d -= min(d);

    // calculate the wavy film height

    const scalar mean = mean_->value(t);
    const scalar a = a_->value(t);
    const scalar omega = omega_->value(t);

    operator==(mean + a*sin(omega*constant::mathematical::twoPi*d));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::wavyFilmHeightFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    mean_->writeData(os);
    a_->writeData(os);
    omega_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wavyFilmHeightFvPatchScalarField
    );
}


// ************************************************************************* //
