/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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

#include "flowRateInletVelocityFuelConvertCoupledFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
//#include "regionProperties.H"
#include "basicThermo.H"
#include "surfaceFields.H"

#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::
flowRateInletVelocityFuelConvertCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nbrPhiName_("none"),
    phiName_("phi"),
    rhoName_("rho"),
    hocPyr_(0.0)
{}


Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::
flowRateInletVelocityFuelConvertCoupledFvPatchVectorField
(
    const flowRateInletVelocityFuelConvertCoupledFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    hocPyr_(ptf.hocPyr_)
{}


Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::
flowRateInletVelocityFuelConvertCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    nbrPhiName_(dict.lookupOrDefault<word>("nbrPhi", "phi")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    hocPyr_(readScalar(dict.lookup("hocPyr")))
{}


Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::
flowRateInletVelocityFuelConvertCoupledFvPatchVectorField
(
    const flowRateInletVelocityFuelConvertCoupledFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    hocPyr_(ptf.hocPyr_)
{}


Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::
flowRateInletVelocityFuelConvertCoupledFvPatchVectorField
(
    const flowRateInletVelocityFuelConvertCoupledFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    hocPyr_(ptf.hocPyr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
//    const mappedPatchBase& distMap = mpp;

    scalarList phi =
        nbrPatch.lookupPatchField<surfaceScalarField, scalar>(nbrPhiName_);


    // get heat of combustion of the gasious fuel
    const basicThermo& thermo =
        db().lookupObject<basicThermo>("thermophysicalProperties");

    const singleStepReactingMixture<gasHThermoPhysics>& singleMixture
    (
        dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&>
        //refCast<const singleStepReactingMixture<gasThermoPhysics> >
        (thermo)
    );

    // heat of combustion [J/kg]
    scalar qFuel(singleMixture.qFuel().value());



    scalarField hocPyrData(nbrPatch.size(),hocPyr_);



    // convert to equivalent gaseous fuel
    //phi = phi * hocPyr_ / qFuel;
    phi = phi * hocPyrData / qFuel;
   
    mpp.distribute(phi);

    const surfaceScalarField& phiName =
        db().lookupObject<surfaceScalarField>(phiName_);


    // a simpler way of doing this would be nice
    //scalar avgU = -flowRate_/gSum(patch().magSf());
    scalarField U(-phi/patch().magSf());

    vectorField n(patch().nf());

//    const surfaceScalarField& phi =
//        db().lookupObject<surfaceScalarField>(phiName_);

    if (phiName.dimensions() == dimVelocity*dimArea)
    {
        // volumetric flow-rate
        operator==(n*U);
    }
    else if (phiName.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // mass flow-rate
        operator==(n*U/rhop);

        if (debug)
        {
            scalar phi = gSum(rhop*(*this) & patch().Sf());
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrMesh.name() << ':'
                << nbrPatch.name() << ':'
                << this->internalField().name() << " :"
                << " mass flux[Kg/s]:" << -phi
                << endl;
        }
    }
    else
    {
        FatalErrorIn
        (
            "flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::flowRateInletVelocityFuelConvertCoupledFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeKeyword("nbrPhi") << nbrPhiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("hocPyr") << hocPyr_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletVelocityFuelConvertCoupledFvPatchVectorField
   );
}


// ************************************************************************* //
