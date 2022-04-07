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

#include "filmHeightInletVelocityWavyFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Random.H" //kvm
#include "constants.H" //kvm

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmHeightInletVelocityWavyFvPatchVectorField::
filmHeightInletVelocityWavyFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    deltafName_("deltaf")
{}


Foam::filmHeightInletVelocityWavyFvPatchVectorField::
filmHeightInletVelocityWavyFvPatchVectorField
(
    const filmHeightInletVelocityWavyFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    deltafName_(ptf.deltafName_)
{}


Foam::filmHeightInletVelocityWavyFvPatchVectorField::
filmHeightInletVelocityWavyFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    deltafName_(dict.lookupOrDefault<word>("deltaf", "deltaf"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::filmHeightInletVelocityWavyFvPatchVectorField::
filmHeightInletVelocityWavyFvPatchVectorField
(
    const filmHeightInletVelocityWavyFvPatchVectorField& fhivpvf
)
:
    fixedValueFvPatchVectorField(fhivpvf),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


Foam::filmHeightInletVelocityWavyFvPatchVectorField::
filmHeightInletVelocityWavyFvPatchVectorField
(
    const filmHeightInletVelocityWavyFvPatchVectorField& fhivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fhivpvf, iF),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityWavyFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const fvPatchField<scalar>& deltafp =
        patch().lookupPatchField<volScalarField, scalar>(deltafName_);

    vectorField n(patch().nf());

    static label seed=0;
    static Random ranGen_(seed);//kvm
    const scalar twoPi=2.0*constant::mathematical::pi;
    scalar rand=0.0;
    scalar time=this->db().time().timeOutputValue();
    const pointField& centers = this->patch().patch().faceCentres();
    scalar spatialFrequency_=32.0;
    scalar amplitude_=1.0;
    forAll(n,i){
        ranGen_.randomise01(rand); 
        rand-=0.5; 
        n[i][1]=0.0+amplitude_*sin(spatialFrequency_*twoPi*centers[i][1]);
    }
    /*Info << n ;*/
    const scalarField& magSf = patch().magSf();

    operator==(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + ROOTVSMALL));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmHeightInletVelocityWavyFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("deltaf", "deltaf", deltafName_);
    writeEntry("value",os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityWavyFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        filmHeightInletVelocityWavyFvPatchVectorField
    );
}


// ************************************************************************* //
