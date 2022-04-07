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

#include "filmHeightInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "surfaceFilmModel.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
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


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& ptf,
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


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    deltafName_(dict.lookupOrDefault<word>("deltaf", "deltaf"))
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf
)
:
    fixedValueFvPatchVectorField(fhivpvf),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fhivpvf, iF),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    // delatf has not yet been solved for, 
    // therefore initial deltafp value will be "value" from deltaf bc file
    const fvPatchField<scalar>& deltafp =
        patch().lookupPatchField<volScalarField, scalar>(deltafName_);

    vectorField n(patch().nf());
    const scalarField& magSf = patch().magSf();

    static bool first=true;
    if(!first){
        first=false;
//###################Nusselt velocity######################
        if(1){
            typedef const Foam::regionModels::surfaceFilmModels::kinematicSingleLayer filmModelType;

            filmModelType& film =
                refCast<filmModelType>(db().time().lookupObject
                                       <Foam::regionModels::surfaceFilmModels::surfaceFilmRegionModel>
                                       ("surfaceFilmProperties")
                    );
            volVectorField gTan2(film.gTan());
            const label patchI = patch().index();
            gTan2.boundaryField()[patchI];
            const fvPatchField<scalar>& mup =
                patch().lookupPatchField<volScalarField, scalar>("muf");
            scalarField nu(deltafp.size()); 
            nu=mup/rhop;
            scalarField gT(deltafp.size()); 
            gT=mag(gTan2.boundaryField()[patchI]);
            /*operator==(0.1*n);*/
            /*Info << "gT " << gT<<endl;*/
            Info << "deltafp " << deltafp<<endl;
            Info << "deltafp.size " << deltafp.size() << endl;
            /*Info << "nu " << nu<<endl;*/
            /*Info << "n " << n<<endl;*/
            operator==(-(gT*deltafp*deltafp/3.0/nu)*n);
        }
//###################Nusselt velocity######################

    }
    else{
        // DEBUG(deltafp[0]);
        // DEBUG(phip[0]);
        // DEBUG(rhop[0]);
        // DEBUG(magSf[0]);
        // const vectorField u(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + ROOTVSMALL));
        // DEBUG(u[0][1]);
        operator==(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + ROOTVSMALL));
        /*operator==(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + SMALL));*/
    }
    

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmHeightInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("deltaf", "deltaf", deltafName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::operator=
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
        filmHeightInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
