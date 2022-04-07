/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "inclinedFilmNusseltHeightFvPatchScalarField.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    GammaMean_(),
    a_(),
    omega_()
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const inclinedFilmNusseltHeightFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    GammaMean_(ptf.GammaMean_.clone()),
    a_(ptf.a_.clone()),
    omega_(ptf.omega_.clone())
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    filmRegionName_
    (
        dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    GammaMean_(Function1<scalar>::New("GammaMean", dict)),
    a_(Function1<scalar>::New("a", dict)),
    omega_(Function1<scalar>::New("omega", dict))
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const inclinedFilmNusseltHeightFvPatchScalarField& wmfrhpsf
)
:
    fixedValueFvPatchScalarField(wmfrhpsf),
    filmRegionName_(wmfrhpsf.filmRegionName_),
    GammaMean_(wmfrhpsf.GammaMean_.clone()),
    a_(wmfrhpsf.a_.clone()),
    omega_(wmfrhpsf.omega_.clone())
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const inclinedFilmNusseltHeightFvPatchScalarField& wmfrhpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wmfrhpsf, iF),
    filmRegionName_(wmfrhpsf.filmRegionName_),
    GammaMean_(wmfrhpsf.GammaMean_.clone()),
    a_(wmfrhpsf.a_.clone()),
    omega_(wmfrhpsf.omega_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inclinedFilmNusseltHeightFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve the film region from the database

    const regionModels::regionModel& region =
        db().time().lookupObject<regionModels::regionModel>(filmRegionName_);

    const regionModels::surfaceFilmModels::kinematicSingleLayer& film =
        dynamic_cast
        <
            const regionModels::surfaceFilmModels::kinematicSingleLayer&
        >(region);

    // Calculate the vector tangential to the patch

    // Note: normal pointing into the domain
    const vectorField n(-patch().nf());

    const volVectorField& nHat = film.nHat();

    const vectorField nHatp(nHat.boundaryField()[patchi].patchInternalField());

    vectorField nTan(nHatp ^ n);
    nTan /= mag(nTan) + ROOTVSMALL;

    // Calculate distance in patch tangential direction

    const vectorField& Cf = patch().Cf();
    scalarField d(nTan & Cf);

    // Calculate the wavy film height

    const scalar t = db().time().timeOutputValue();

    // TODO: currently re-evaluating the entire gTan field to return this patch
    const scalarField gTan(film.gTan()().boundaryField()[patchi] & n);

    if (patch().size() && (max(mag(gTan)) < SMALL))
    {
        WarningIn
        (
            "void Foam::inclinedFilmNusseltHeightFvPatchScalarField::"
            "updateCoeffs()"
        )
            << "Tangential gravity component is zero.  This boundary condition "
            << "is designed to operate on patches inclined with respect to "
            << "gravity"
            << endl;
    }

    const volScalarField& mu = film.mu();
    // The use of the first internal cell for properties doesn't seem correct, kvm
    const scalarField mup(mu.boundaryField()[patchi].patchInternalField());
    // const scalarField mup(mu.boundaryField()[patchi]);

    const volScalarField& rho = film.rho();
    const scalarField rhop(rho.boundaryField()[patchi].patchInternalField());
    // const scalarField rhop(rho.boundaryField()[patchi]);
    
    // calculate the wavy film height

    const scalar GMean = GammaMean_->value(t);
    const scalar a = a_->value(t);
    const scalar omega = omega_->value(t);


    // solve for deltaMean via Bi-Section method
    scalar fxC = 10.0;
    scalar tol = 0.00001;
    label iter=0;
    scalar deltaMeanA = 2e-2;
    scalar deltaMeanB = 2e-6;
    scalar deltaMeanC;

    //  Info << " acoeff " << a << nl;
    scalarField C(pow(3.0*sqr(mup/rhop)/mup/(gTan + ROOTVSMALL), 0.33333333));
    scalarField delta(deltaMeanA + a*sin(omega*constant::mathematical::twoPi*d));

    while(fabs(fxC)>tol){

        label deltaSize = delta.size();
        reduce(deltaSize, sumOp<label>());
        delta = (deltaMeanA + a*sin(omega*constant::mathematical::twoPi*d));
        scalar fxA = GMean - gSum(pow(delta/C,3))/deltaSize; 
        // DEBUG(delta.size());
        // DEBUG(deltaSize);
        if(fxA > 0.0){
            WarningIn
                (
                 "void Foam::inclinedFilmNusseltHeightFvPatchScalarField::"
                 "updateCoeffs()"
                )
                << "Initial guess for deltaMeanA too low:"
                << deltaMeanA
                << nl;
        }


        delta = (deltaMeanB + a*sin(omega*constant::mathematical::twoPi*d));
        scalar fxB = GMean - gSum(pow(delta/C,3))/deltaSize;

        if(fxB < 0.0){
            WarningIn
                (
                 "void Foam::inclinedFilmNusseltHeightFvPatchScalarField::"
                 "updateCoeffs()"
                )
                << "Initial guess for deltaMeanB too high:"
                << deltaMeanB
                << nl;
        }
        //  Info << "fxA " << fxA << nl;
        //  Info << "fxB " << fxB << nl;

        deltaMeanC = 0.5*(deltaMeanA+deltaMeanB);

        delta = (deltaMeanC + a*sin(omega*constant::mathematical::twoPi*d));
        fxC = GMean - gSum(pow(delta/C,3))/deltaSize;

        //  Info << iter << " fxC " << fxC << nl;

        if( fxC < 0.0 ) {
            deltaMeanA = deltaMeanC;
        }
        else{
            deltaMeanB = deltaMeanC;
        }
        iter++;
        if(iter>1e5){
            WarningIn
                (
                 "void Foam::inclinedFilmNusseltHeightFvPatchScalarField::"
                 "updateCoeffs()"
                )
                << "Maximum number of bisection method iterations reached."
                << nl;
            break;
        }
    }

    scalar deltaMean = deltaMeanC;
    delta = (deltaMean + a*sin(omega*constant::mathematical::twoPi*d));
    //  Info << "delta " << delta << nl;

    // scalarField G(GMean + GMean*a*sin(omega*constant::mathematical::twoPi*d));
    
    // correction to mean 
    // scalar uncorrectedMean = gSum(G)/G.size();
    // DEBUG(uncorrectedMean);
    // G *= GMean/(uncorrectedMean+SMALL);

    // set internal values
    // volScalarField& rho2 = const_cast<volScalarField&>(film.rho());
    // const labelUList& faceCells = patch().faceCells();
    // forAll(faceCells,i)
    // {
    //     label cellI = faceCells[i];
    //     rho2[cellI] = 1.0;
    // }
    // Info << rho2 << nl;

    scalarField G(pow(delta/C,3));

    const scalarField Re(max(G, scalar(0.0))/mup);

    // DEBUG(mup[0]);
    // DEBUG(rhop[0]);
    // DEBUG(gTan[0]);
    // DEBUG(Re[0]);
    // DEBUG(G[0]);
    // const scalarField h(pow(3.0*sqr(mup/rhop)/(-gTan + ROOTVSMALL), 0.33333333)*pow(Re, 0.33333333));
    // DEBUG(h[0]);

    operator==
    (
        cbrt(3*sqr(mup/rhop)/(gTan + ROOTVSMALL))*cbrt(Re)
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


//void Foam::inclinedFilmNusseltHeightFvPatchScalarField::updateCoeffs2()
//{
//    if (updated())
//    {
//        return;
//    }
//
//    const label patchi = patch().index();
//
//    const scalar t = db().time().timeOutputValue();
//
//    // retrieve the film region from the database
//
//    const regionModels::regionModel& region =
//        db().time().lookupObject<regionModels::regionModel>
//        (
//            "surfaceFilmProperties"
//        );
//
//    const regionModels::surfaceFilmModels::kinematicSingleLayer& film =
//        dynamic_cast
//        <
//            const regionModels::surfaceFilmModels::kinematicSingleLayer&
//        >(region);
//
//    // calculate the vector tangential to the patch
//
//    const vectorField n(patch().nf());
//
//    const volVectorField& nHat = film.nHat();
//
//    const vectorField nHatp(nHat.boundaryField()[patchi].patchInternalField());
//
//    vectorField nTan(nHatp ^ n);
//    nTan /= mag(nTan) + ROOTVSMALL;
//
//    // calculate distance in patch tangential direction
//
//    const vectorField& Cf = patch().Cf();
//    scalarField d(nTan & Cf);
//
//    // calculate the wavy film height
//
//    const scalar GMean = GammaMean_->value(t);
//    const scalar a = a_->value(t);
//    const scalar omega = omega_->value(t);
//
//    scalarField G(GMean + GMean*a*sin(omega*constant::mathematical::twoPi*d));
//    
//    // correction to mean 
//    scalar uncorrectedMean = gSum(G)/G.size();
//    // DEBUG(uncorrectedMean);
//    G *= GMean/(uncorrectedMean+SMALL);
//
//    const volScalarField& mu = film.mu();
//    // The use of the first internal cell for properties doesn't seem correct, kvm
//    const scalarField mup(mu.boundaryField()[patchi].patchInternalField());
//    // const scalarField mup(mu.boundaryField()[patchi]);
//
//    const volScalarField& rho = film.rho();
//    const scalarField rhop(rho.boundaryField()[patchi].patchInternalField());
//    // const scalarField rhop(rho.boundaryField()[patchi]);
//    
//    // set internal values
//    // volScalarField& rho2 = const_cast<volScalarField&>(film.rho());
//    // const labelUList& faceCells = patch().faceCells();
//    // forAll(faceCells,i)
//    // {
//    //     label cellI = faceCells[i];
//    //     rho2[cellI] = 1.0;
//    // }
//    // Info << rho2 << nl;
//
//    const scalarField Re(max(G, 0.0)/mup);
//
//    // TODO: currently re-evaluating the entire gTan field to return this patch
//    const scalarField gTan(film.gTan()().boundaryField()[patchi] & n);
//
//    if (patch().size() && (max(mag(gTan)) < SMALL))
//    {
//       WarningIn
//        (
//            "void Foam::inclinedFilmNusseltHeightFvPatchScalarField::"
//            "updateCoeffs()"
//        )
//            << "Tangential gravity component is zero.  This boundary condition "
//            << "is designed to operate on patches inclined with respect to "
//            << "gravity"
//            << nl;
//    }
//
//    // DEBUG(mup[0]);
//    // DEBUG(rhop[0]);
//    // DEBUG(gTan[0]);
//    // DEBUG(Re[0]);
//    // DEBUG(G[0]);
//    // const scalarField h(pow(3.0*sqr(mup/rhop)/(-gTan + ROOTVSMALL), 0.33333333)*pow(Re, 0.33333333));
//    // DEBUG(h[0]);
//
//    operator==
//    (
//        pow(3.0*sqr(mup/rhop)/(-gTan + ROOTVSMALL), 0.33333333)*pow(Re, 0.33333333)
//    );
//
//    fixedValueFvPatchScalarField::updateCoeffs();
//}
//
//
void Foam::inclinedFilmNusseltHeightFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    GammaMean_->writeData(os);
    a_->writeData(os);
    omega_->writeData(os);
    writeEntry("value2", os); //kvm, value written twice if this line active
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inclinedFilmNusseltHeightFvPatchScalarField
    );
}


// ************************************************************************* //
