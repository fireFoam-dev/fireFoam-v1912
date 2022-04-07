/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenCFD Ltd.
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

#include "corrugatedCardboardFreeburnFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "turbulenceModel.H"
//#include "thermoSingleLayer.H"
#include "pyrolysisModel.H"

#include "constants.H"

#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//const corrugatedCardboardFreeburnFvPatchScalarField::filmModelType&
//corrugatedCardboardFreeburnFvPatchScalarField::
//filmModel() const
//{
//    HashTable<const filmModelType*> models
//        = db().time().lookupClass<filmModelType>();
//
//    forAllConstIter(HashTable<const filmModelType*>, models, iter)
//    {
//        if (iter()->regionMesh().name() == filmRegionName_)
//        {
//            return *iter();
//        }
//    }
//
//
//    FatalErrorIn
//    (
//        "const corrugatedCardboardFreeburnFvPatchScalarField::"
//        "filmModelType& "
//        "corrugatedCardboardFreeburnFvPatchScalarField::"
//        "filmModel() const"
//    )
//        << "Unable to locate film region " << filmRegionName_
//        << abort(FatalError);
//
//    return **models.begin();
//}


const corrugatedCardboardFreeburnFvPatchScalarField::
pyrolysisModelType&
corrugatedCardboardFreeburnFvPatchScalarField::
pyrModel() const
{
    HashTable<const pyrolysisModelType*> models =
        db().time().lookupClass<pyrolysisModelType>();

    forAllConstIter(HashTable<const pyrolysisModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == pyrolysisRegionName_)
        {
            return *iter();
        }
    }

    FatalErrorIn
    (
        "const corrugatedCardboardFreeburnFvPatchScalarField::"
        "pyrolysisModelType& "
        "corrugatedCardboardFreeburnFvPatchScalarField::"
        "pyrModel() const"
    )
        << "Unable to locate pyrolysis region " << pyrolysisRegionName_
        << abort(FatalError);

    return **models.begin();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

corrugatedCardboardFreeburnFvPatchScalarField::
corrugatedCardboardFreeburnFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase
    (
        patch(),
        "undefined",
        "undefined",
        "undefined-K",
        "undefined-alpha"
    ),
    //filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("undefined-Tnbr"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    neighbourFieldConvectiveName_("undefined-neigbourFieldConvectiveName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    fieldConvectiveName_("undefined-fieldConvectiveName"),
    KName_("undefined-K"),
    convectiveScaling_(1.0),
    convectiveCoefficient_(1.0),
    //filmDeltaDry_(0.0),
    //filmDeltaWet_(0.0),
    qBound_(false),
    qMax_(2.0e5),
    //qExtra_(0),
    oldMode_(unknown)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


corrugatedCardboardFreeburnFvPatchScalarField::
corrugatedCardboardFreeburnFvPatchScalarField
(
    const corrugatedCardboardFreeburnFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    //filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    fieldConvectiveName_(psf.fieldConvectiveName_),
    KName_(psf.KName_),
    convectiveScaling_(psf.convectiveScaling_),
    convectiveCoefficient_(psf.convectiveCoefficient_),
    //filmDeltaDry_(psf.filmDeltaDry_),
    //filmDeltaWet_(psf.filmDeltaWet_),
    qBound_(psf.qBound_),
    qMax_(psf.qMax_),
    //qExtra_(psf.qExtra_),
    oldMode_(psf.oldMode_)
{}


corrugatedCardboardFreeburnFvPatchScalarField::
corrugatedCardboardFreeburnFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    //filmRegionName_
    //(
    //    dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    //),
    pyrolysisRegionName_
    (
        dict.lookupOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookup("Tnbr")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    neighbourFieldConvectiveName_(dict.lookup("neighbourFieldConvectiveName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    fieldConvectiveName_(dict.lookup("fieldConvectiveName")),
    KName_(dict.lookup("K")),
    convectiveScaling_(dict.lookupOrDefault<scalar>("convectiveScaling", 1.0)),
    convectiveCoefficient_(dict.lookupOrDefault<scalar>("convectiveCoefficient", 1.0)),
    //filmDeltaDry_
    //(
    //        dict.lookupOrDefault<scalar>("filmDeltaDry", 0.0000)
    //),
    //filmDeltaWet_
    //(
    //        dict.lookupOrDefault<scalar>("filmDeltaWet", 0.0002)
    //),
    qBound_(dict.lookupOrDefault<bool>("qBound", false)),
    qMax_(dict.lookupOrDefault<scalar>("qMax", 2.0e5)),
    //qExtra_(dict.lookupOrDefault<scalar>("qExtra", 0.0)),

    oldMode_(unknown)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "corrugatedCardboardFreeburnFvPatchScalarField::"
            "corrugatedCardboardFreeburnFvPatchScalarField\n"
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


corrugatedCardboardFreeburnFvPatchScalarField::
corrugatedCardboardFreeburnFvPatchScalarField
(
    const corrugatedCardboardFreeburnFvPatchScalarField&
        psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    //filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    fieldConvectiveName_(psf.fieldConvectiveName_),
    KName_(psf.KName_),
    convectiveScaling_(psf.convectiveScaling_),
    convectiveCoefficient_(psf.convectiveCoefficient_),
    //filmDeltaDry_(psf.filmDeltaDry_),
    //filmDeltaWet_(psf.filmDeltaWet_),
    qBound_(psf.qBound_),
    qMax_(psf.qMax_),
    //qExtra_(psf.qExtra_),
    oldMode_(psf.oldMode_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void corrugatedCardboardFreeburnFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchI = patch().index();
    const label nbrPatchI = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    scalarField intFld(patchInternalField());

    const corrugatedCardboardFreeburnFvPatchScalarField&
        nbrField =
        refCast
        <
            const corrugatedCardboardFreeburnFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField& Tp = *this;

    const scalarField K(this->kappa(*this));
    const scalarField nbrK(nbrField.kappa(*this));

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);

    scalarField nbrConvFlux(KDeltaNbr*(intFld - nbrIntFld));

    scalarField nbrTotalFlux = nbrConvFlux;
    // scalarList nbrRadField(nbrPatch.size(), 0.0);
    //scalarList QrCoupled(nbrPatch.size(), 0.0);//kvm
    //scalarList Tfilm(nbrPatch.size(), 0.0);//kvm
    //scalarList alpha(nbrPatch.size(), 0.0);//kvm
    //scalarList filmConv(nbrPatch.size(), 0.0);//kvm
    //scalarList filmDelta(nbrPatch.size(), 0.0);//kvm
    //scalarList myRadField(patch().size(), 0.0);
    //scalarList qExtraZone(nbrPatch.size(), 0.0);//NR
    

    scalarList radiField(nbrPatch.size(), 0.0);//Ning
    scalarList convField(nbrPatch.size(), 0.0);//Ning
    //scalarList radiField(nbrPatch.size(), 0.0);//Ning

    const pyrolysisModelType& pyrolysis = pyrModel();
    //const filmModelType& film = filmModel();
    scalarField Twall(patch().size(),0.0);

    // In solid
    if(neighbourFieldRadiativeName_ != "none") //nbr Radiation Qr
    {

        scalarField nbrConvFlux(convectiveCoefficient_*(intFld - nbrIntFld));

        //const label filmPatchI =
        //    pyrolysis.nbrCoupledPatchID(film, patchI);

        //const scalarField Qconvw(film.qconvw(filmPatchI));

        //// kvm, Qconvw is not right
        //filmConv =
        //    pyrolysis.mapRegionPatchField<scalar>
        //    (
        //        film,
        //        "qFilmToWall",
        //        patchI,
        //        true
        //    );

        //QrCoupled =
        //    pyrolysis.mapRegionPatchField<scalar>
        //    (
        //        film,
        //        "qin",
        //        patchI,
        //        true
        //    );

        //Tfilm =
        //    pyrolysis.mapRegionPatchField<scalar>
        //    (
        //        film,
        //        "Tf",
        //        patchI,
        //        true
        //    );

        //filmDelta =
        //    pyrolysis.mapRegionPatchField<scalar>
        //    (
        //        film,
        //        "deltaf",
        //        patchI,
        //        true
        //    );
        //alpha =
        //    pyrolysis.mapRegionPatchField<scalar>
        //    (
        //        film,
        //        "alpha",
        //        patchI,
        //        true
        //    );
        //bool qExtraZoneExists = false;
        //if(film.regionMesh().foundObject<volScalarField>("correctionZone"))
        //{
        //    qExtraZone =
        //    pyrolysis.mapRegionPatchField<scalar>
        //    (
        //        film,
        //        "correctionZone",
        //        patchI,
        //        true
        //    );
        //    Info<<"Add extra heat flux at MLR front: "<<qExtra_<<endl;
        //    qExtraZoneExists = true;
        //}

            const fvMesh& mesh = patch().boundaryMesh().mesh();
            if (! (mesh.foundObject<radiation::radiationModel>("radiationProperties")))
            {
                FatalErrorIn
                (
                    "corrugatedCardboardFreeburnFvPatchScalarField::"
                    "corrugatedCardboardFreeburnFvPatchScalarField\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const DimensionedField<scalar, volMesh>& iF,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "\n    radiationProperties file not found in pyrolysis region\n"
                    << exit(FatalError);
            }
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

            scalarField tabsorptivity
            (
                radiation.absorptionEmission().a()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );


        // nbrRadField =
        //     nbrPatch.lookupPatchField<volScalarField, scalar>
        //     (
        //         neighbourFieldRadiativeName_
        //     );

        //convField =	//Ning
        //    nbrPatch.lookupPatchField<volScalarField, scalar>
        //    (
        //        neighbourFieldConvectiveName_
        //    );
        radiField =	//Ning
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                "qin"
            );

	    mpp.distribute(radiField);	//Ning
        convField =	//Ning
            nbrPatch.lookupPatchField<surfaceScalarField, scalar>
            (
                "convectiveHeatFlux_T"
            );

	    mpp.distribute(convField);	//Ning

        //radiField =	//Ning
        //    nbrPatch.lookupPatchField<volScalarField, scalar>
        //    (
        //        neighbourFieldRadiativeName_
        //    );

	    //mpp.distribute(radiField);	//Ning

        // Note: the Qr radiative flux is positive outgoing.
        // For a hot solid radiating into a cold fluid Qr will be negative.


        // Swap to obtain full local values of neighbour radiative heat flux
        // field
        // mpp.distribute(nbrRadField);

    //Qrad is negative going out of the solid
    //Qcon is positive going out of the solid

        if (debug)
        {
            // scalar Qr = gSum(nbrRadField*patch().magSf());

            Info<< mesh.name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " :" << nl
                // << "     radiative heat  [W] : " << Qr << nl
//kvm                << "     predicted wallT [K] : " << gAverage(Twall) << nl
                << endl;
        }

        label nFixed = 0;

         //// Estimate wetness of the film (1: wet , 0: dry)
         //scalarField ratio
         //(
         //   min
         //   (
         //       max
         //       (
         //           (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
         //           scalar(0.0)
         //       ),
         //       scalar(1.0)
         //   )
         //);

	//tmp<scalarField> kappaTmp = (kappa(*this));
	//fvPatchField kappaTmp(kappa(*this));

        forAll(*this, i)
        {
            //scalar qConvWet = -filmConv[i];
            ////scalar qConvDry = nbrConvFlux[i];
            //scalar qConvDry = -convField[i]; //Ning

            const scalar sigma = constant::physicoChemical::sigma.value();
            //scalar qRadWet = 0.0; // all film absorption takes place in film model
            //// use of 'ratio' leads to non-conservation of energy
            //// if use 'ratio', then need to use it in film model absorption as well
            ////kvm scalar qRadWet =-QrCoupled[i];
            //scalar qRadDry =
            //   //-temissivity[i]*QrCoupled[i]
            //   -tabsorptivity[i]*QrCoupled[i]
            //   //-tabsorptivity[i]*radiField[i] //Ning
            //   //kvm -(1.0-ratio)*temissivity[i]*QrCoupled[i]
            //   +temissivity[i]*sigma*pow4(operator[](i));

            //scalar qConv = alpha[i]*qConvWet + (1.0-alpha[i])*qConvDry;
            //scalar qRad  = alpha[i]*qRadWet  + (1.0-alpha[i])*qRadDry;

            ////nbrTotalFlux[i] = qConv + qRad;
            //scalar extraHeatFlux(0);
            //if(qExtraZoneExists && qExtraZone[i] > 0.5)
            //{
            //    //extraHeatFlux = -qExtra_;
            //    extraHeatFlux = (alpha[i]-1.0)*qExtra_;
            //    //Info<<"facei :"<<i<<tab<<qConv<<tab<<qRad<<tab<<extraHeatFlux<<endl;
            //}
            scalar qConv = -convField[i];
            scalar qRad = -tabsorptivity[i]*radiField[i] + temissivity[i]*sigma*pow4(operator[](i));
            nbrTotalFlux[i] = qConv + qRad;

	        if(qBound_)
	        {
		        if(nbrTotalFlux[i] > qMax_)
		        {
		            nbrTotalFlux[i] = qMax_;
		        }
		        else if(nbrTotalFlux[i] < -qMax_)
		        {
		            nbrTotalFlux[i] = -qMax_;
		        }
	        }
            this->refValue()[i] = operator[](i);  // not used
            //this->refGrad()[i] = -nbrTotalFlux[i]/kappa(*this)()[i];
            this->refGrad()[i] = -nbrTotalFlux[i]/K[i];
            //DEBUG(nbrTotalFlux[i]);
            //DEBUG(kappa(*this)()[i]);

            this->valueFraction()[i] = 0.0;
            nFixed++;
        }

        if (debug)
        {
            Pout<< "Using " << nFixed << " fixedValue out of " << this->size()
                << endl;
        }
    }
    else // In fluid
    {
        Twall = nbrIntFld;

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        //scalar Qc = gSum(nbrConvFlux*patch().magSf());
        scalar Qc = gSum(convField*patch().magSf()); //Ning
        // scalar Qr = gSum(radField*patch().magSf());
        scalar Qt = gSum(nbrTotalFlux*patch().magSf());

        // Info<< mesh.name() << ':'
        //     << patch().name() << ':'
        //     << this->internalField().name() << " -> "
        //     << nbrMesh.name() << ':'
        //     << nbrPatch.name() << ':'
        //     << this->internalField().name() << " :"
        //     << " heatFlux:" << Qc
        //     // << " radiativeFlux:" << Qr
        //     << " totalFlux:" << Qt
        //     << " walltemperature "
        //     << " min:" << gMin(*this)
        //     << " max:" << gMax(*this)
        //     << " avg:" << gAverage(*this)
        //     << endl;
    }
}


void corrugatedCardboardFreeburnFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    //os.writeEntryIfDifferent<word>
    //(
    //    "filmRegion",
    //    "surfaceFilmProperties",
    //    filmRegionName_
    //);
    os.writeEntryIfDifferent<word>
    (
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldConvectiveName")<< //Ning
        neighbourFieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldConvectiveName")<< //Ning
        fieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<< //Ning
        fieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("K")<< //Ning
        KName_ << token::END_STATEMENT << nl;
    //os.writeKeyword("qExtra")<< //Ning
    //    qExtra_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    corrugatedCardboardFreeburnFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
