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

#include "thermoSingleLayerExtraQ.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"

#include "pyrolysisModel.H"

#include "fvc.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayerExtraQ, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerExtraQ, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void thermoSingleLayerExtraQ::updateSubmodels()
{
    thermoSingleLayerFmBox::updateSubmodels();

    getHeatFluxCorrectionZone();
}

void thermoSingleLayerExtraQ::getHeatFluxCorrectionZone()
{
    phiSolid_.correctBoundaryConditions();
    qConvect_.correctBoundaryConditions();
    qRadR_.correctBoundaryConditions();

    totalHeatFlux_ = qConvect_+qRadR_;
    massLossRate_ = phiSolid_;

    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        const label patchi = intCoupledPatchIDs_[i];
        scalarField& phiSolidW = phiSolid_.boundaryFieldRef()[patchi];
        scalarField& qcW  = qConvect_.boundaryFieldRef()[patchi];
        scalarField& qrW  = qRadR_.boundaryFieldRef()[patchi];
        //const fvPatch& patch = regionMesh().boundary()[patchi];
        const labelUList& faceCells = regionMesh().boundary()[patchi].faceCells();
        forAll(phiSolidW,facei)
        {
            totalHeatFlux_[faceCells[facei]] = qcW[facei]+qrW[facei];
            massLossRate_[faceCells[facei]] = phiSolidW[facei];
        }
    }

    massLossRate_.correctBoundaryConditions();

    //forAll(totalHeatFlux_,cellI)
    //{
    //    totalHeatFlux_[cellI] = qConvect_[cellI]+qRadR_[cellI];
    //}
    //totalHeatFlux_.correctBoundaryConditions();

    const vectorField& Cv = regionMesh().C();
    const scalarField& cellVolume = regionMesh().V();
    scalar mlrCritical(1);

    forAll(massLossRate_, cellI)
    {
        correctionZone_[cellI] = 0;

        if(massLossRate_[cellI] > 3.0*mlrCritical)
        {
            if(burntZone_[cellI] == 0)
            {
                burntZone_[cellI] = 1.0;
            }
        }

        if(massLossRate_[cellI] < mlrCritical)
        {
            const labelList& nbrCL = regionMesh().cellCells()[cellI];
            forAll(nbrCL,nI)
            {
                label nbrCellID = nbrCL[nI];
                if((massLossRate_[nbrCellID] > mlrCritical) && (burntZone_[cellI] == 0))
                {
                    correctionZone_[cellI] = 1;
                }
            }
        }
    }

    forAll(massLossRate_.boundaryField(),pI)
    {
        const fvPatchField<scalar>& mlrB = massLossRate_.boundaryField()[pI];
        if(mlrB.type() == "processor")
        {
            const labelUList& faceCells = mlrB.patch().faceCells();
            forAll(mlrB,fi)
            {
                label cellO = faceCells[fi];
                if((burntZone_[cellO] == 0) && (massLossRate_[cellO] < mlrCritical) && (mlrB[fi] > mlrCritical))
                {
                    correctionZone_[cellO] = 1;
                }
            }
        }
    }

    correctionZone_.correctBoundaryConditions();
    burntZone_.correctBoundaryConditions();

    for(int ipc = 1; ipc < preHeatingCells_; ipc++)
    {
        Info<<"CorrectionZone No.: "<<ipc+1<<endl;
        forAll(massLossRate_, cellI)
        {
            if((burntZone_[cellI] == 0) && (correctionZone_[cellI] == 0) && (massLossRate_[cellI] < mlrCritical))
            {
                const labelList& nbrCL = regionMesh().cellCells()[cellI];
                forAll(nbrCL,nI)
                {
                    label nbrCellID = nbrCL[nI];
                    if(correctionZone_[nbrCellID] == ipc)
                    {
                        correctionZone_[cellI] = 1 + ipc;
                        //Info<<"on-cell: "<<Cv[cellI].x()<<tab<<Cv[cellI].y()<<tab<<Cv[cellI].z()<<endl;
                    }
                }
            }
        }
        //forAll(correctionZone_.boundaryField(),pI)
        //{
        //    const fvPatchField<scalar>& correctionZoneB = correctionZone_.boundaryField()[pI];
        //    if(correctionZoneB.type() == "processor")
        //    {
        //        const labelUList& faceCells = correctionZoneB.patch().faceCells();
        //        forAll(correctionZoneB,fi)
        //        {
        //            label cellO = faceCells[fi];
        //            if((burntZone_[cellO] == 0) && (correctionZone_[cellO] == 0) && (massLossRate_[cellO] < mlrCritical))
        //            {
        //                if(correctionZone_[fi] == ipc)
        //                {
        //                    correctionZone_[cellO] = 1 + ipc;
        //                    //Info<<"on-face: "<<Cv[cellO].x()<<tab<<Cv[cellO].y()<<tab<<Cv[cellO].z()<<endl;
        //                }
        //            }
        //        }
        //    }
        //}
        correctionZone_.correctBoundaryConditions();
    }

    forAll(correctionZone_, cellI)
    {
        if(Cv[cellI].z() > preHeatingZoneMaxHeight_)
        {
            correctionZone_[cellI] = 0;
        }
    }
    correctionZone_.correctBoundaryConditions();

    const scalar deltaT = time_.deltaTValue();
    forAll(delta_, cellI)
    {
        scalar dWater(0.0);
        if(delta_[cellI] > 0.00005)
        {
            dWater = absorptionRate_*deltaT;
        }
        else
        {
            if((waterAbsorbed_[cellI] > 0) && (delta_[cellI] < 0.00002))
            {
                dWater = min(0.0,-totalHeatFlux_[cellI]*deltaT/2.5e6);
            }
        }
        waterAbsorbed_[cellI] = min(saturationThickness_, waterAbsorbed_[cellI]+dWater);
    }
    waterAbsorbed_.correctBoundaryConditions();

        //scalarField& phiSolidW = phiSolid_.boundaryFieldRef()[patchi];
        //scalarField& qcW  = qConvect_.boundaryFieldRef()[patchi];
        //scalarField& qrW  = qRadR_.boundaryFieldRef()[patchi];
        //scalarField& totalHeatFluxW  = totalHeatFlux_.boundaryFieldRef()[patchi];
        //scalarField& massLossRateW  = massLossRate_.boundaryFieldRef()[patchi];
        ////scalarField& qrW  = qRadR_.boundaryFieldRef()[patchi];

        //const fvPatch& patch = regionMesh().boundary()[patchi];
        //const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        //(
        //    patch.patch()
        //);
        //const polyMesh& nbrMesh = mpp.sampleMesh();
        //const fvPatch& nbrPatch = refCast<const fvMesh>
        //(
        //    nbrMesh
        //).boundary()[mpp.samplePolyPatch().index()];

        //fvPatchScalarField phiW = nbrPatch.lookupPatchField<surfaceScalarField, scalar>("phi");
        //fvPatchScalarField qcTW = nbrPatch.lookupPatchField<surfaceScalarField, scalar>("convectiveHeatFlux_T");
        //mpp.distribute(phiW);
        //mpp.distribute(qcTW);
        //phiSolidW = -phiW;
        //qcW = qcTW;
        //massLossRateW = phiSolidW/regionMesh().boundary()[patchi].magSf();
        //totalHeatFluxW = qcW + qrW;

        //forAll(phiW, facei)
        //{
        //    massLossRate_[facei] = -phiW[cellI]/regionMesh().boundary()[patchI].magSf()[cellI];
        //    totalHeatFlux_[facei] = qRadR_[facei]+qcW[facei];
        //}


        //label patchI = intCoupledPatchIDs_[i];
        //label patchIp = primaryPatchIDs_[i];
        //const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
        //scalarList flowRate(pp.faceCells().size(), 0.0);
        ////const scalarField& phiw = primaryMesh().phi().boundaryField()[patchI]; //boundary()[patchI].phi();
        //const scalarField& phiw = primaryMesh().lookupObject<surfaceScalarField>("phi").boundaryField()[patchIp];
        //const scalarField& qcw = primaryMesh().lookupObject<surfaceScalarField>("convectiveHeatFlux_T").boundaryField()[patchIp];
        //Info<<primaryPatchIDs_<<tab<<intCoupledPatchIDs_<<tab<<phiw.size()<<tab<<flowRate.size()<<endl;
        //const scalarField& phiw = primaryMesh().patch().lookupPatchField<surfaceScalarField, scalar>("phi");
        //const regionModels::regionModel& gasRegion =
        //        db().time().lookupObject<regionModels::regionModel>
        //    (
        //        "region0"
        //    );

        //typedef regionModels::pyrolysisModels::pyrolysisModel
        //        pyrolysisModelType;
        //const regionModels::regionModel& pyrolysisRegion =
        //        db().time().lookupObject<regionModels::regionModel>
        //    (
        //        "pyrolysisProperties"
        //    );
        //const pyrolysisModelType& pyrolysisModel =
        //        dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
        //pyrolysisModelType& pyrolysis =
        //    const_cast<pyrolysisModelType&>(pyrolysisModel);

        //scalarList MLR(pp.faceCells().size(), 0.0);
        //MLR = mapRegionPatchInternalField<scalar>
        //        (
        //            pyrolysis,
        //            "phiGas",
        //            patchI,
        //            true
        //        );

        //scalarList tempList(pp.faceCells().size(), 0.0);
        //tempList = mapRegionPatchInternalField<scalar>
        //        (
        //            pyrolysis,
        //            "T",
        //            patchI,
        //            true
        //        );

        //scalarList tempMLR(pp.faceCells().size(), 0.0);
        //tempMLR = mapRegionPatchField<scalar>
        //        (
        //            //pyrolysis,
        //            //primaryMesh(),
        //            gasRegion,
        //            "phi",
        //            patchI,
        //            true
        //        );

        //forAll(flowRate, cellI)
        //{
        //    massLossRate_[cellI] = -phiw[cellI]/regionMesh().boundary()[patchI].magSf()[cellI];
        //    totalHeatFlux_[cellI] = qcw[cellI] + qRadR_.boundaryField()[patchI][cellI];
        //    //Info<<"dbg-- "<<cellI<<phiw[cellI];
        //    //massLossRate_[cellI] = MLR[cellI]/regionMesh().boundary()[patchI].magSf()[cellI];
        //    //flowRate[cellI] = -phiSolid_.boundaryField()[patchI][cellI];
        //    //massLossRate_[cellI] = flowRate[cellI]/regionMesh().boundary()[patchI].magSf()[cellI];
        //    //totalHeatFlux_[cellI] = qConvect_.boundaryField()[patchI][cellI]+qRadR_.boundaryField()[patchI][cellI];
        //    //TsolidSurface_[cellI] = tempList[cellI];
        //}
        
    //}

    //massLossRate_.correctBoundaryConditions();

    ////- Count the page difference at the pyrolysis front
    //forAll(filmPaperID_, cellI)
    //{
    //    filmPositiveDelNeighbor_[cellI] = 0;
    //    filmNegtiveDelNeighbor_[cellI]  = 0;
    //    const labelList& nbrCL = regionMesh().cellCells()[cellI];
    //    forAll(nbrCL,nI)
    //    {
    //        label nbrCellID = nbrCL[nI];
    //        scalar deltP = filmPaperID_[nbrCellID] - filmPaperID_[cellI];
    //        if(deltP > 0)
    //        {
    //            filmPositiveDelNeighbor_[cellI]++;
    //            vector vAverage = 0.5*(Cv[cellI] + Cv[nbrCellID]);
    //            scalar dAverage = mag(Cv[cellI] - vAverage);
    //            if(dAverage < filmMinPyrolysisDistance_[cellI])
    //            {
    //                filmMinPyrolysisDistance_[cellI] = dAverage;
    //                filmPyrolysisFrontCellPosition_[cellI] = vAverage;
    //                filmPyrolysisFrontPaperID_[cellI] = filmPaperID_[cellI];
    //            }
    //        }
    //        else if(deltP < 0)
    //        {
    //            filmNegtiveDelNeighbor_[cellI]++;
    //        }
    //        else
    //        {
    //            //- Do nothing
    //        }
    //    }
    //}

    ////- Continue counting (for parallel run), update processor boundary
    //forAll(filmPaperID_.boundaryField(),pI)
    //{
    //    const fvPatchField<scalar>& paper_f = filmPaperID_.boundaryField()[pI];
    //if(paper_f.type() == "processor")
    //{
    //    //Info<<"patchType-processor: "<<paper_f.type()<<endl;
    //        const labelUList& faceCells = paper_f.patch().faceCells();
    //        vectorField faceB = paper_f.patch().Cf();
    //        forAll(paper_f, fi)
    //        {
    //            label cellO = faceCells[fi];
    //        scalar deltP = paper_f[fi] - filmPaperID_[cellO];
    //        if(deltP > 0)
    //        {
    //            filmPositiveDelNeighbor_[cellO]++;
    //            scalar dnbr = mag(Cv[cellO] - faceB[fi]);
    //            if(dnbr < filmMinPyrolysisDistance_[cellO])
    //            {
    //                filmMinPyrolysisDistance_[cellO] = dnbr;
    //                filmPyrolysisFrontCellPosition_[cellO] = faceB[fi];
    //                filmPyrolysisFrontPaperID_[cellO] = filmPaperID_[cellO];
    //            }
    //        }
    //        else if(deltP < 0)
    //        {
    //            filmNegtiveDelNeighbor_[cellO]++;
    //        }
    //        else
    //        {
    //            //- Do nothing
    //        }
    //        }
    //}
    //}
    //forAll(filmPaperID_, cellI)
    //{
    //if(filmNegtiveDelNeighbor_[cellI] > 0)
    //{
    //    filmMinPyrolysisDistance_[cellI] = 1000;
    //    filmPyrolysisFrontCellPosition_[cellI] = vector(0,0,-1000.0);
    //    filmPyrolysisFrontPaperID_[cellI] = -filmPaperID_[cellI];
    //}
    //}
    //
    //// - update nearest pyrolysis front distance from neighbour cell
    //for(int ipy=0;ipy<5;ipy++)
    //{
    //    forAll(filmPyrolysisFrontCellPosition_,cellI)
    //    {
    //        //- For cells away from the pyrolysis front
    //    if(
    //    (filmPositiveDelNeighbor_[cellI] == 0) 
    //    && 
    //    (filmNegtiveDelNeighbor_[cellI] == 0)
    //      )
    //        {
    //            const labelList& nbrCL = regionMesh().cellCells()[cellI];
    //            forAll(nbrCL, nI)
    //            {
    //                label nbrCellID = nbrCL[nI];
    //                if(filmPaperID_[cellI] == filmPyrolysisFrontPaperID_[nbrCellID])
    //                {
    //                    vector vnbr = Cv[cellI] - filmPyrolysisFrontCellPosition_[nbrCellID];
    //                    scalar dist = mag(vnbr);
    //                    if(dist < filmMinPyrolysisDistance_[cellI])
    //                    {
    //                        filmMinPyrolysisDistance_[cellI] = dist;
    //                        filmPyrolysisFrontCellPosition_[cellI] = 
    //                filmPyrolysisFrontCellPosition_[nbrCellID];
    //                        filmPyrolysisFrontPaperID_[cellI] = filmPaperID_[cellI];
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    //- For parallel run, update processor boundary
    //    filmPyrolysisFrontCellPosition_.correctBoundaryConditions();
    //    filmPyrolysisFrontPaperID_.correctBoundaryConditions();

    //    forAll(filmPyrolysisFrontCellPosition_.boundaryField(),pI)
    //    {
    //        const fvPatchField<vector>& pyCpmin_f = 
    //        filmPyrolysisFrontCellPosition_.boundaryField()[pI];
    ////if(pyCpmin_f.type() == "processor")
    ////{
    //        const fvPatchField<scalar>& pyPage_f = 
    //        filmPyrolysisFrontPaperID_.boundaryField()[pI];
    //        const labelUList& faceCells = pyCpmin_f.patch().faceCells();
    //        vectorField faceB = pyCpmin_f.patch().Cf();
    //        forAll(pyCpmin_f, fi)
    //        {
    //            label cellO = faceCells[fi];
    //            //- For cells away from the pyrolysis front
    //        if(
    //            (filmPositiveDelNeighbor_[cellO] == 0) 
    //            && 
    //            (filmNegtiveDelNeighbor_[cellO] == 0)
    //          )
    //            {
    //                vector vnbr = Cv[cellO] - pyCpmin_f[fi];
    //                if(filmPaperID_[cellO] == pyPage_f[fi])
    //                {
    //                    scalar dist = mag(vnbr);
    //                    if(dist < filmMinPyrolysisDistance_[cellO])
    //                    {
    //                        filmMinPyrolysisDistance_[cellO] = dist;
    //                        filmPyrolysisFrontCellPosition_[cellO] = pyCpmin_f[fi];
    //                        filmPyrolysisFrontPaperID_[cellO] = filmPaperID_[cellO];
    //                    }
    //                }
    //            }
    //        }
    //    //}
    //}
    //}

    ////- Correct bad cells 
    //forAll(filmPaperID_, cellI)
    //{
    //if(fabs(filmPyrolysisFrontPaperID_[cellI]) != filmPaperID_[cellI])
    //    {
    //    Info<<"Correct-BAD-cells: "<<cellI<<tab<<filmPaperID_[cellI]<<tab
    //        <<filmPyrolysisFrontPaperID_[cellI]<<tab
    //    <<filmMinPyrolysisDistance_[cellI]<<tab
    //    <<filmPyrolysisFrontCellPosition_[cellI]<<endl;
    //        filmMinPyrolysisDistance_[cellI] = 1000;
    //        filmPyrolysisFrontCellPosition_[cellI] = vector(0,0,-1000.0);
    //    //filmPyrolysisFrontPaperID_[cellI] = filmPaperID_[cellI];
    //    //pyPage_[cellI] = paper_[cellI];
    //    }
    //}

    ////- Set thermally-thin zone based on minimum distance to pyFront
    //forAll(filmMinPyrolysisDistance_,cellI)
    //{
    //    scalar halfCellGridSize = 0.5*sqrt(cellVolume[cellI]*1000.0);
    //scalar delPy = filmPeelingDistance_ - filmMinPyrolysisDistance_[cellI];
    //if(delPy + halfCellGridSize < 0)
    //{
    //    filmPaperPeelingZone_[cellI] = 0.0;
    //}
    //else if(delPy < halfCellGridSize)
    //{
    //    filmPaperPeelingZone_[cellI] = (1.0 + delPy/halfCellGridSize)/2.0;
    //}
    //else
    //{
    //    filmPaperPeelingZone_[cellI] = 1.0;
    //}
    //}

    ////- Special treatment when pyrolysis fronts are stacked
    //forAll(filmPaperID_, cellI)
    //{
    //    if(
    //        (filmPositiveDelNeighbor_[cellI] > 0) 
    //        && 
    //        (filmNegtiveDelNeighbor_[cellI] > 0)
    //      )
    //    {
    //        filmPaperPeelingZone_[cellI] = 1.0;
    //    }

	//    //- Turn on paper stack effect when paper sheet reach a threshold
	//    if(
	//        (filmPaperID_[cellI] > paperStackEffectID_) 
	//        && 
	//        (filmNegtiveDelNeighbor_[cellI] > 0)
	//      )
	//    {
	//        filmPaperPeelingZone_[cellI] = 0.0;
	//    }
    //}

    //// - Calculate maximum thermally thin region height
    //const vectorField& cellCentre = regionMesh().cellCentres();
    //filmPaperPyrolysisFront_ = dimensionedScalar("zero",dimLength,0.0);
    //forAll(filmPaperPeelingZone_, cellI)
    //{
    //    //scalar cellGridSize = sqrt(cellVolume[cellI]*1000.0);
    //    if(filmPaperPeelingZone_[cellI] > 0.5)
    //    {
    //        filmPaperPyrolysisFront_[cellI] = cellCentre[cellI].z();
    //    }
    //}

    //scalar maxPyHeight = max(filmPaperPyrolysisFront_).value();
    //Info<<"ThinPaper-MaximumHeight: "<<time_.value()<<tab<<maxPyHeight<<endl;

    //if(maxPyHeight < filmMinPyrolysisHeight_)
    //{
    //    Info<<"Turnoff thermally thin model!"<<endl;
    //    forAll(filmPaperPeelingZone_, cellI)
    //    {
    //        filmPaperPeelingZone_[cellI] = 0;
    //    }
    //}

    //filmTpaper_.correctBoundaryConditions();
    //filmTsurface_.correctBoundaryConditions();
    //filmPaperPeelingZone_.correctBoundaryConditions();
    //filmQnet_.correctBoundaryConditions();
    //filmBlockFactor_.correctBoundaryConditions();
    //filmMinPyrolysisDistance_.correctBoundaryConditions();
    //filmMassReleaseRate_.correctBoundaryConditions();
    //filmPositiveDelNeighbor_.correctBoundaryConditions();
    //filmNegtiveDelNeighbor_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerExtraQ::thermoSingleLayerExtraQ
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    //thermoSingleLayerFmBox(modelType, mesh, g, regionType, readFields)
    thermoSingleLayerFmBox(modelType, mesh, g, regionType, readFields),
    //preHeatingDistance_(coeffs_.lookupOrDefault<scalar>("preHeatingDistance", 0.0254)),
    preHeatingCells_(coeffs_.lookupOrDefault<scalar>("preHeatingCells", 1.0)),
    preHeatingZoneMaxHeight_(coeffs_.lookupOrDefault<scalar>("preHeatingZoneMaxHeight", 1.5)),
    absorptionRate_(coeffs_.lookupOrDefault<scalar>("absorptionRate", 0.005)),
    saturationThickness_(coeffs_.lookupOrDefault<scalar>("SaturationThickness", 1.0)),
    massLossRate_
    (
        IOobject
        (
            "solidMLR",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    totalHeatFlux_
    (
        IOobject
        (
            "Qtotal",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    correctionZone_
    (
        IOobject
        (
            "correctionZone",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    burntZone_
    (
        IOobject
        (
            "burntZone",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            //IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    waterAbsorbed_
    (
        IOobject
        (
            "waterAbsorbed",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            //IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    phiSolid_
    (
        IOobject
        (
            "mlrSolid", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        //dimensionedScalar("zero", dimMass/dimTime, 0.0)
        dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    qConvect_
    (
        IOobject
        (
            //"convectiveHeatFlux_T", // Same name as T on primary region to enable mapping
            "QcWallFunction",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        //dimensionedScalar("zero", dimensionSet(1, 0, -3, 0 , 0, 0, 0), 0.0)
        dimensionedScalar("zero", dimensionSet(1, 0, -3, 0 , 0, 0, 0), 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    qRadR_
    (
        IOobject
        (
            "qr", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimensionSet(1, 0, -3, 0 , 0, 0, 0), 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    )
{
    Info<<"Add extra heat flux"<<endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerExtraQ::~thermoSingleLayerExtraQ()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
