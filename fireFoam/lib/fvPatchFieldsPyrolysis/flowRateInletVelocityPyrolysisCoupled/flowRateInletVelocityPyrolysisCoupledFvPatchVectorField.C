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

#include "flowRateInletVelocityPyrolysisCoupledFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "basicThermo.H"
#include "surfaceFields.H"

#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"

#include "reactingMixture.H"
#include "constIsoSolidTransport.H"
#include "hConstThermo.H"
#include "rhoConst.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::
flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nbrPhiName_("none"),
    phiName_("phi"),
    rhoName_("rho"),
    hocSolid_(0.0)
{}


Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::
flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
(
    const flowRateInletVelocityPyrolysisCoupledFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    hocSolid_(ptf.hocSolid_)
{}


Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::
flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
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
    hocSolid_(readScalar(dict.lookup("hocSolid")))
{}


Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::
flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
(
    const flowRateInletVelocityPyrolysisCoupledFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    hocSolid_(ptf.hocSolid_)
{}


Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::
flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
(
    const flowRateInletVelocityPyrolysisCoupledFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    hocSolid_(ptf.hocSolid_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::updateCoeffs()
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

    scalarList phi =
        nbrPatch.lookupPatchField<surfaceScalarField, scalar>(nbrPhiName_);

    // kvm, for debugging restart
    scalar phiMax = -VGREAT;
    scalar phiMin =  VGREAT;
    phiMax = max(phi);
    phiMin = min(phi);
    reduce(phiMax,maxOp<scalar>());
    reduce(phiMin,minOp<scalar>());
    Info << "flowRateInletVelocityPyrolysisCoupled::phiMin: " << tab << db().time().timeName() << tab << phiMin << endl;
    Info << "flowRateInletVelocityPyrolysisCoupled::phiMax: " << tab << db().time().timeName() << tab << phiMax << endl;

    // get heat of combustion of the gasious fuel
    const basicThermo& thermo =
        db().lookupObject<basicThermo>("thermophysicalProperties");

    const singleStepReactingMixture<gasHThermoPhysics>& singleMixture
    (
        dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&>
        (thermo)
    );

    // heat of combustion [J/kg]
    scalar qFuel(singleMixture.qFuel().value());

    // get solidThermo from mapped pyrolysis region
    //const basicSolidThermo& solidThermo =
    //    mpp.sampleMesh().lookupObject<basicSolidThermo>("solidThermophysicalProperties");
    
    const basicThermo& solidThermo =
        mpp.sampleMesh().lookupObject<basicThermo>("thermophysicalProperties");

    // access density of char and v
    const reactingMixture<constIsoSolidTransport<species::thermo<hConstThermo<rhoConst<specie> >, sensibleEnthalpy > > >& mcSolidMixture
    (
        dynamic_cast<const reactingMixture<constIsoSolidTransport<species::thermo<hConstThermo<rhoConst<specie> >, sensibleEnthalpy > > >&>
        (solidThermo)
    );

    const label charIndex = mcSolidMixture.species()["char"];
    const label vIndex = mcSolidMixture.species()["v"];
    //const label vIndex = mcSolidMixture.species()["wood"];

    const scalar rhoChar = mcSolidMixture.rho(charIndex,1,298.0);  //plug in any P and T for rhoConst
    const scalar rhoV = mcSolidMixture.rho(vIndex,1,298.0); 

    //Info << "density of char is " << rhoChar << endl;
    //Info << "density of v is " << rhoV << endl;

    // Heat of combustion of char
    scalar hocChar = 32.8e6; //[W/kg]
        
    // Heat of combustion of gaseous pyrolysate 
    //scalar hocPyr = hocSolid_ * rhoV / (rhoV - rhoChar);    //no Char Oxi
    scalar hocPyr = (hocSolid_ * rhoV - hocChar * rhoChar) / (rhoV - rhoChar);
    //Info << "HOC for pyrolysate is " << hocPyr << endl;

    scalarField hocPyrData(nbrPatch.size(),hocPyr);


    // convert to equivalent gaseous fuel
    // scalar phiMax = max(max(phi),-VGREAT);
    // reduce(phiMax,maxOp<scalar>());
    // Info << "flowRateInletVelocityPyrolysisCoupled::phiMax: " << tab << db().time().timeName() << tab << phiMax << endl;
    //phi = phi * hocPyr / qFuel;
    phi = phi * hocPyrData / qFuel;

    mpp.distribute(phi);

    const surfaceScalarField& phiName =
        db().lookupObject<surfaceScalarField>(phiName_);

    scalarField U(-phi/patch().magSf());

    vectorField n(patch().nf());

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
            scalar phi(gSum(rhop*(*this) & patch().Sf()));
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
            "flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::flowRateInletVelocityPyrolysisCoupledFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeKeyword("nbrPhi") << nbrPhiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("hocSolid") << hocSolid_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletVelocityPyrolysisCoupledFvPatchVectorField
   );
}


// ************************************************************************* //
