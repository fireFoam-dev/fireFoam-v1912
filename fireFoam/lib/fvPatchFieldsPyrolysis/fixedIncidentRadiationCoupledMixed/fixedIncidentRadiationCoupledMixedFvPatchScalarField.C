/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "fixedIncidentRadiationCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
//#include "regionProperties.H"
#include "basicThermo.H"
#include "LESModel.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedIncidentRadiationCoupledMixedFvPatchScalarField::
fixedIncidentRadiationCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBaseFF(p, "undefined", "undefined", scalarField::null(), scalarField::null()),
    neighbourFieldName_("undefined-neighbourFieldName"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    KName_("undefined-K"),
    QrIncident_(0.0)
//    emissivity_(p.size(), 0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


fixedIncidentRadiationCoupledMixedFvPatchScalarField::
fixedIncidentRadiationCoupledMixedFvPatchScalarField
(
    const fixedIncidentRadiationCoupledMixedFvPatchScalarField& ptf,
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
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourFieldRadiativeName_(ptf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(ptf.fieldRadiativeName_),
    KName_(ptf.KName_),
    QrIncident_(ptf.QrIncident_)
//    emissivity_(ptf.emissivity_)
{}


fixedIncidentRadiationCoupledMixedFvPatchScalarField::
fixedIncidentRadiationCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBaseFF(p, dict),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    KName_(dict.lookup("K")),
    QrIncident_(readScalar(dict.lookup("QrIncident")))
//    emissivity_(p.size(), 0.0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "fixedIncidentRadiationCoupledMixedFvPatchScalarField::"
            "fixedIncidentRadiationCoupledMixedFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


fixedIncidentRadiationCoupledMixedFvPatchScalarField::
fixedIncidentRadiationCoupledMixedFvPatchScalarField
(
    const fixedIncidentRadiationCoupledMixedFvPatchScalarField&
        wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    radiationCoupledBaseFF
    (
        wtcsf.patch(),
        wtcsf.emissivityMethod(),
        wtcsf.absorptivityMethod(),
        wtcsf.emissivity_,
        wtcsf.absorptivity_
    ),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    neighbourFieldRadiativeName_(wtcsf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(wtcsf.fieldRadiativeName_),
    KName_(wtcsf.KName_),
    QrIncident_(wtcsf.QrIncident_)
//    emissivity_(wtcsf.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
fixedIncidentRadiationCoupledMixedFvPatchScalarField::K() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (KName_ == "none")
    {
        const compressible::turbulenceModel& turbModel =
            db().lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    internalField().group()
                )
            );

        tmp<volScalarField> talpha = turbModel.alphaEff();

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        return
            talpha().boundaryField()[patch().index()]
           *thermo.Cp()().boundaryField()[patch().index()];
    }

    else if (KName_ == "solidThermo")
    {
        const fvMesh& mesh = patch().boundaryMesh().mesh();
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>("thermophysicalProperties");

        //scalarField K_ = thermo.kappa(patch().index());
        return thermo.kappa(patch().index());
    }

    else if (mesh.objectRegistry::foundObject<volScalarField>(KName_))
    {
        return patch().lookupPatchField<volScalarField, scalar>(KName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KName_))
    {
        const symmTensorField& KWall =
            patch().lookupPatchField<volSymmTensorField, scalar>(KName_);

        vectorField n(patch().nf());

        return n & KWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "fixedIncidentRadiationCoupledMixedFvPatchScalarField"
            "::K() const"
        )   << "Did not find field " << KName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'K' to 'none', 'solidThermo', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void fixedIncidentRadiationCoupledMixedFvPatchScalarField::
updateCoeffs()
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
    //const mapDistribute& distMap = mpp.map();
    //const mappedPatchBase& distMap = mpp;

    scalarField intFld(patchInternalField());

    const fixedIncidentRadiationCoupledMixedFvPatchScalarField&
    nbrField =  refCast
        <const fixedIncidentRadiationCoupledMixedFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldName_
            )
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarList radField(nbrPatch.size(),0.0);  
    scalarField Twall(patch().size(),0.0);

    // In solid
    if(neighbourFieldRadiativeName_ != "none") //nbr Radiation Qr
    {
        radField =
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldRadiativeName_
            );
   
        // Swap to obtain full local values of neighbour radiative heat flux field
        mpp.distribute(radField);

//        emissivity_ =
//            patch().lookupPatchField<volScalarField, scalar>("emissivity");

//        const scalarField temissivity = emissivity();
            const fvMesh& mesh = patch().boundaryMesh().mesh();
            const radiation::radiationModel& radiation =
                mesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField temissivity
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );

        scalarField nbrTotalFlux(-radField);



/*
        Twall =
            (radField + myKDelta*intFld + nbrKDelta*nbrIntFld)
            /(myKDelta + nbrKDelta);
*/

        //hard code for now.
        const scalar sigma_ = 5.67e-8;

        scalarField gradT_((QrIncident_*temissivity-temissivity*sigma_*pow(intFld,4))/K());

        forAll(*this, i)
        {
		//fixed gradient BC, use internal T to replace Tb for re-radiation
                this->refValue()[i] = operator[](i);  // not used
                this->refGrad()[i] = gradT_[i];
                this->valueFraction()[i] = 0.0;
        }

    }
    else // In fluid
    {
        radField =
            patch().lookupPatchField<volScalarField, scalar>
            (
                fieldRadiativeName_
            );

        Twall = nbrIntFld;

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qr = gSum(radField*patch().magSf());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " radiativeFlux:" << Qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void fixedIncidentRadiationCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<< fieldRadiativeName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("QrIncident")<< QrIncident_
        << token::END_STATEMENT << nl;
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedIncidentRadiationCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
