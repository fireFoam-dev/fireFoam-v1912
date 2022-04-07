/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "massFlowInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Random.H" //kvm
#include "mathematicalConstants.H" //kvm

#include "surfaceFilmModel.H"
#include "kinematicSingleLayer.H"

#include <iostream>
#include <fstream>
using namespace std;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massFlowInletFvPatchScalarField::
massFlowInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    phiName_("phi"),
    rhoName_("rhof"),
    muName_("muf"),
    deltafName_("deltaf"),
    gamma_(0.1)
{}


Foam::massFlowInletFvPatchScalarField::
massFlowInletFvPatchScalarField
(
    const massFlowInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    muName_(ptf.muName_),
    deltafName_(ptf.deltafName_),
    gamma_(ptf.gamma_)
{}


Foam::massFlowInletFvPatchScalarField::
massFlowInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    filmRegionName_
    (
        dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rhof")),
    muName_(dict.lookupOrDefault<word>("mu", "muf")),
    deltafName_(dict.lookupOrDefault<word>("deltaf", "deltaf")),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 0.1))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::massFlowInletFvPatchScalarField::
massFlowInletFvPatchScalarField
(
    const massFlowInletFvPatchScalarField& fhivpvf
)
:
    fixedValueFvPatchScalarField(fhivpvf),
    filmRegionName_(fhivpvf.filmRegionName_),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    muName_(fhivpvf.muName_),
    deltafName_(fhivpvf.deltafName_),
    gamma_(fhivpvf.gamma_)
{}


Foam::massFlowInletFvPatchScalarField::
massFlowInletFvPatchScalarField
(
    const massFlowInletFvPatchScalarField& fhivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fhivpvf, iF),
    filmRegionName_(fhivpvf.filmRegionName_),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    muName_(fhivpvf.muName_),
    deltafName_(fhivpvf.deltafName_),
    gamma_(fhivpvf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massFlowInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const fvPatchField<scalar>& mup =
        patch().lookupPatchField<volScalarField, scalar>(muName_);

    const fvPatchField<scalar>& deltafp =
        patch().lookupPatchField<volScalarField, scalar>(deltafName_);

    scalarField dfa(deltafp.size()); 
    scalarField Re(deltafp.size()); 
    scalarField nu(deltafp.size()); 
    scalarField filmThickness(deltafp.size()); 
    scalarField gT(deltafp.size()); 

    // Retrieve the film region from the database
    const regionModels::regionModel& region =
        db().time().lookupObject<regionModels::regionModel>(filmRegionName_);
    const regionModels::surfaceFilmModels::kinematicSingleLayer& film =
        dynamic_cast
        <
            const regionModels::surfaceFilmModels::kinematicSingleLayer&
        >(region);

    const label patchI = patch().index();
    volVectorField gTan2(film.gTan());
    gTan2.boundaryField()[patchI];

    /*Nusselt solution*/
    Re=gamma_/mup;
    nu=mup/rhop;
    gT=mag(gTan2.boundaryField()[patchI]);
    filmThickness=pow(3*pow(nu,2)/gT,0.3333)*pow(Re,0.3333);
    operator==(filmThickness);

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::massFlowInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("mu", "mu", muName_);
    os.writeEntryIfDifferent<word>("deltaf", "deltaf", deltafName_);
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::massFlowInletFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& pvf
)
{
    fvPatchField<scalar>::operator=(pvf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        massFlowInletFvPatchScalarField
    );
}


// ************************************************************************* //
