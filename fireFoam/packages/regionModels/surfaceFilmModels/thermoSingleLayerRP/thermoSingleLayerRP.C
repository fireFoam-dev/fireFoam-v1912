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

#include "thermoSingleLayerRP.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"

// Sub-models
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "massAbsorptionModel.H"
#include "filmRadiationModel.H"
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

defineTypeNameAndDebug(thermoSingleLayerRP, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerRP, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



void thermoSingleLayerRP::updateSurfaceTemperatures()
{
    thermoSingleLayer::updateSurfaceTemperatures();

    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchI = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchI];

        //-Ning: get roll paper page number
        if(rollPaperModel_)
        {
            Info<<"Using Roll-Paper Pyrolysis Model!"<<endl;
                typedef regionModels::pyrolysisModels::pyrolysisModel
                    pyrolysisModelType;
                const regionModels::regionModel& pyrolysisRegion =
                    db().time().lookupObject<regionModels::regionModel>
                        (
                            "pyrolysisProperties"
                        );
                const pyrolysisModelType& pyrolysisModel =
                    dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
                pyrolysisModelType& pyrolysis =
                    const_cast<pyrolysisModelType&>(pyrolysisModel);

            scalarList pageList(pp.faceCells().size(), 0.0);
            pageList = mapRegionPatchInternalField<scalar>
                    (
                        pyrolysis,
                        "page",
                        patchI,
                        true
                    );
            scalarList tempList(pp.faceCells().size(), 0.0);
            tempList = mapRegionPatchInternalField<scalar>
                    (
                        pyrolysis,
                        "T",
                        patchI,
                        true
                    );
            scalarList tempSurfList(pp.faceCells().size(), 0.0);
            tempSurfList = mapRegionPatchInternalField<scalar>
                    (
                        pyrolysis,
                        "Tsurface",
                        patchI,
                        true
                    );
            scalarList tempQ(pp.faceCells().size(), 0.0);
            tempQ = mapRegionPatchInternalField<scalar>
                    (
                        pyrolysis,
                        "Qnet",
                        patchI,
                        true
                    );
            scalarList tempBlock(pp.faceCells().size(), 0.0);
            tempBlock = mapRegionPatchInternalField<scalar>
                    (
                        pyrolysis,
                        "blockFactor",
                        patchI,
                        true
                    );
            //Info<<"patch size: "<<pageList.size()<<tab<<paper_.size()<<endl;
            forAll(paper_, cellI)
            {
                paper_[cellI] = pageList[cellI];
                Tpaper_[cellI] = tempList[cellI];
                TpaperSurf_[cellI] = tempSurfList[cellI];
                pflux_[cellI] = tempQ[cellI];
                pblock_[cellI] = tempBlock[cellI];
            }
        }
    }

    //-Ning: update boundary condition
    if(rollPaperModel_)
    {
    paper_.correctBoundaryConditions();
    const vectorField& Cv = regionMesh().C();
    const scalarField& cellVolume = regionMesh().V();
    forAll(paper_, cellI)
    {
        const labelList& nbrCL = regionMesh().cellCells()[cellI];
        pyState_[cellI] = 0;
        forAll(nbrCL,nI)
        {
            if(pyState_[cellI] != -1.0)        //No negtive neighbor
            {
            label nbrCellID = nbrCL[nI];
                scalar deltP = paper_[nbrCellID] - paper_[cellI];
                if(deltP > 0)        //Neighbor delaminated
                {
                    pyState_[cellI]++;
            vector vAverage = 0.5*(Cv[cellI] + Cv[nbrCellID]);
            scalar dAverage = mag(Cv[cellI] - vAverage);
            if(dAverage < pyRpmin_[cellI])
            {
                pyRpmin_[cellI] = dAverage;
                pyCpmin_[cellI] = vAverage;
                pyPage_[cellI] = paper_[cellI];
            }
                }
                else if(deltP < 0)        //Cell delaminated, neighbor not
                {
                    pyState_[cellI] = -1.0;
            pyRpmin_[cellI] = 1000;
            pyCpmin_[cellI] = vector(0,0,-1000.0);
            pyPage_[cellI] = -paper_[cellI];
                }
            }
        }
    }
    forAll(paper_.boundaryField(),pI)
    {
        const fvPatchField<scalar>& paper_f = paper_.boundaryField()[pI];
        const labelUList& faceCells = paper_f.patch().faceCells();
        vectorField faceB = paper_f.patch().Cf();
        forAll(paper_f, fi)
        {
        label cellO = faceCells[fi];
        if(pyState_[cellO] != -1.0)           //No negtive neighbor
        {
            scalar deltP = paper_f[fi] - paper_[cellO];
            if(deltP > 0)        //Neighbor delaminated
                {
                    pyState_[cellO]++;
            scalar dnbr = mag(Cv[cellO] - faceB[fi]);
            if(dnbr < pyRpmin_[cellO])
            {
                pyRpmin_[cellO] = dnbr;
                pyCpmin_[cellO] = faceB[fi];
                pyPage_[cellO] = paper_[cellO];
            }
                }
                else if(deltP < 0)        //Cell delaminated, neighbor not
                {
                    pyState_[cellO] = -1.0;
            pyRpmin_[cellO] = 1000;
            pyCpmin_[cellO] = vector(0,0,-1000.0);
            pyPage_[cellO] = -paper_[cellO];
                }
        }
        }
    }
    //pyCpmin_.correctBoundaryConditions();
    //- Possibly another bug
    //-Todo: check the page number of the pyrolysis front
    //-    minimum is asigned only when Pcurrent equal p_pyfront.
    
    // - update nearest pyrolysis front distance from neighbour cell
    for(int ipy=0;ipy<5;ipy++)
    {
        forAll(pyCpmin_,cellI)
        {
            if(pyState_[cellI] == 0)
            {
                const labelList& nbrCL = regionMesh().cellCells()[cellI];
                forAll(nbrCL, nI)
                {
                    label nbrCellID = nbrCL[nI];
                if(paper_[cellI] == pyPage_[nbrCellID])
                {
                    vector vnbr = Cv[cellI] - pyCpmin_[nbrCellID];
                        scalar dist = mag(vnbr);
                if(dist < pyRpmin_[cellI])
                        {
                            pyRpmin_[cellI] = dist;
                            pyCpmin_[cellI] = pyCpmin_[nbrCellID];
                    pyPage_[cellI] = paper_[cellI];
                        }
                }
                }
            }
        }
        pyCpmin_.correctBoundaryConditions();
        pyPage_.correctBoundaryConditions();
        forAll(pyCpmin_.boundaryField(),pI)
        {
            const fvPatchField<vector>& pyCpmin_f = pyCpmin_.boundaryField()[pI];
            const fvPatchField<scalar>& pyPage_f = pyPage_.boundaryField()[pI];
            const labelUList& faceCells = pyCpmin_f.patch().faceCells();
            vectorField faceB = pyCpmin_f.patch().Cf();
            forAll(pyCpmin_f, fi)
            {
            label cellO = faceCells[fi];
                if(pyState_[cellO] == 0)
                {
                    vector vnbr = Cv[cellO] - pyCpmin_f[fi];
                if(paper_[cellO] == pyPage_f[fi])
                {
                scalar dist = mag(vnbr);
                if(dist < pyRpmin_[cellO])
                    {
                        pyRpmin_[cellO] = dist;
                        pyCpmin_[cellO] = pyCpmin_f[fi];
                    pyPage_[cellO] = paper_[cellO];
                    }
                }
                }
            }
        }
    }

    //- Correct bad cells 
    forAll(paper_, cellI)
    {
        if(pyPage_[cellI] < 0)
        {
        pyRpmin_[cellI] = 1000;
        pyCpmin_[cellI] = vector(0,0,-1000.0);
        }
    }
    //forAll(paper_, cellI)
    //{
    //    if((pyState_[cellI] == 0) && (paper_[cellI] != pyPage_[cellI]))
    //    {
    //        Info<<"Bad-Cell: "<<cellI<<endl;
    //        scalar rpMin = 1000.0;
    //    pyRpmin_[cellI] = 0.2;
    //        const labelList& nbrCL = regionMesh().cellCells()[cellI];
    //    forAll(nbrCL, nI)
    //    {
    //        label nbrCellID = nbrCL[nI];
    //        vector vnbr = Cv[cellI] - pyCpmin_[nbrCellID];
    //        scalar dist = mag(vnbr);
    //        if(dist < pyRpmin_[cellI])
    //        {
    //        pyRpmin_[cellI] = dist;
    //        pyCpmin_[cellI] = pyCpmin_[nbrCellID];
    //        }
    //    }
    //    }
    //}

    //- Update pthin based on minimum distance to pyFront
    forAll(pyRpmin_,cellI)
    {
        scalar halfCellGridSize = 0.5*sqrt(cellVolume[cellI]*1000.0);
        vector vCp = Cv[cellI] - pyCpmin_[cellI];
        //- get scaled distance
        scalar scaleVertical;
        if(vCp.z()>0)
        {
        scaleVertical = max(0.1,scaleUp_);
        }
        else
        {
        if(pyState_[cellI] > 1)
        {
            scaleVertical = 1.0;
        }
        else
        {
            scaleVertical = max(0.1,scaleDown_);
        }
        }
        vector vScale = vector(vCp.x(),vCp.y(),vCp.z()/scaleVertical);
        pyRpminScaled_[cellI] = mag(vScale);
        //
        if((pyRpminScaled_[cellI] - halfCellGridSize) < delDistance_)
        {
            if((pyRpminScaled_[cellI] + halfCellGridSize) < delDistance_)
            {
                pthin_[cellI] = 1;
            }
            else
            {
                pthin_[cellI] = 0.5*(
                    (delDistance_ - pyRpminScaled_[cellI])/halfCellGridSize
                    + 1.0);
            }
        }
        else
        {
            pthin_[cellI] = 0;
        }
    }

    const vectorField& cellCentre = regionMesh().cellCentres();
    // - Calculate maximum thermally thin region height
    pthinHeight_ = dimensionedScalar("zero",dimLength,0.0);
    forAll(pthin_, cellI)
    {
        //scalar cellGridSize = sqrt(cellVolume[cellI]*1000.0);
        if(pthin_[cellI] > 0.2)
        {
            pthinHeight_[cellI] = cellCentre[cellI].z();
        }
    }
    scalar maxPyHeight = max(pthinHeight_).value();
    Info<<"ThinPaper-MaximumHeight: "<<time_.value()<<tab<<maxPyHeight<<endl;
    if(maxPyHeight < minPyHeight_)
    {
        Info<<"Turnoff thermally thin model!"<<endl;
        forAll(pthin_, cellI)
        {
            pthin_[cellI] = 0;
        }
    }

    //ps_.correctBoundaryConditions();
    Tpaper_.correctBoundaryConditions();
    TpaperSurf_.correctBoundaryConditions();
    pthin_.correctBoundaryConditions();
    pflux_.correctBoundaryConditions();
    pblock_.correctBoundaryConditions();
    pyState_.correctBoundaryConditions();
    pyRpmin_.correctBoundaryConditions();
    //pyCpmin_.correctBoundaryConditions();
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerRP::thermoSingleLayerRP
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    thermoSingleLayer(modelType, mesh, g, regionType, readFields),
    //rollPaperModel_(readBool(coeffs_.lookupOrDefault("rollPaperModel", false))),
    rollPaperModel_(coeffs_.lookupOrDefault<bool>("rollPaperModel", false)),
    minPyHeight_(coeffs_.lookupOrDefault<scalar>("minPyHeight", 1.0)),
    delDistance_(coeffs_.lookupOrDefault<scalar>("delDistance", 0.025)),
    scaleUp_(coeffs_.lookupOrDefault<scalar>("scaleUp", 1.0)),
    scaleDown_(coeffs_.lookupOrDefault<scalar>("scaleDown", 0.5)),

    //- new model variable to calculate distance
    pyState_
    (
        IOobject
        (
            "pyState",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pyState", dimTime/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    pyRpmin_
    (
        IOobject
        (
            "pyRpmin",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pyRpmin", dimLength, 1000.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    pyRpminScaled_
    (
        IOobject
        (
            "pyRpminScaled",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pyRpminScaled", dimLength, 1000.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    pyCpmin_
    (
        IOobject
        (
            "pyCpmin",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedVector("pyCpmin", dimLength, vector(0,0,-1000)),
        zeroGradientFvPatchScalarField::typeName
    ),

    pyPage_
    (
        IOobject
        (
            "pyPage",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pyPage", dimTime/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    paper_
    (
        IOobject
        (
            "paperNumber",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("paperNumber", dimTime/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    pblock_
    (
        IOobject
        (
            "paperBlock",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("paperBlock", dimTime/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    pflux_
    (
        IOobject
        (
            "netFlux",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("netFlux", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    pthin_
    (
        IOobject
        (
            "pthin",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pthin", dimTime/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    pthinHeight_
    (
        IOobject
        (
            "pthinHeight",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pthinHeight", dimLength, 0.0)
    ),

    Tpaper_
    (
        IOobject
        (
            "paperT",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("paperT", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    TpaperSurf_
    (
        IOobject
        (
            "paperSurfT",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("paperSurfT", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    )

{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerRP::~thermoSingleLayerRP()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
