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

#include "thermoSingleLayerRollPaper.H"
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

defineTypeNameAndDebug(thermoSingleLayerRollPaper, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerRollPaper, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void thermoSingleLayerRollPaper::updateSurfaceTemperatures()
{
    thermoSingleLayer::updateSurfaceTemperatures();

    updateRollPaperThermallyThinZone();
}

void thermoSingleLayerRollPaper::updateRollPaperThermallyThinZone()
{
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchI = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchI];

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
                    "paperID",
                    patchI,
                    true
                );

    //- The following variables are for debug purpose
    //  (tempList, tempSurfList, tempQ, tempBlock, tempMLR)
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

        scalarList tempMLR(pp.faceCells().size(), 0.0);
        tempMLR = mapRegionPatchInternalField<scalar>
                (
                    pyrolysis,
                    "massReleaseRate",
                    patchI,
                    true
                );

        forAll(filmPaperID_, cellI)
        {
            filmPaperID_[cellI] = pageList[cellI];
            filmTpaper_[cellI] = tempList[cellI];
            filmTsurface_[cellI] = tempSurfList[cellI];
            filmQnet_[cellI] = tempQ[cellI];
            filmBlockFactor_[cellI] = tempBlock[cellI];
            filmMassReleaseRate_[cellI] = tempMLR[cellI];
        }
        
    }

    filmPaperID_.correctBoundaryConditions();
    const vectorField& Cv = regionMesh().C();
    const scalarField& cellVolume = regionMesh().V();

    //- Count the page difference at the pyrolysis front
    forAll(filmPaperID_, cellI)
    {
        filmPositiveDelNeighbor_[cellI] = 0;
        filmNegtiveDelNeighbor_[cellI]  = 0;
        const labelList& nbrCL = regionMesh().cellCells()[cellI];
        forAll(nbrCL,nI)
        {
            label nbrCellID = nbrCL[nI];
            scalar deltP = filmPaperID_[nbrCellID] - filmPaperID_[cellI];
            if(deltP > 0)
            {
                filmPositiveDelNeighbor_[cellI]++;
                vector vAverage = 0.5*(Cv[cellI] + Cv[nbrCellID]);
                scalar dAverage = mag(Cv[cellI] - vAverage);
                if(dAverage < filmMinPyrolysisDistance_[cellI])
                {
                    filmMinPyrolysisDistance_[cellI] = dAverage;
                    filmPyrolysisFrontCellPosition_[cellI] = vAverage;
                    filmPyrolysisFrontPaperID_[cellI] = filmPaperID_[cellI];
                }
            }
            else if(deltP < 0)
            {
                filmNegtiveDelNeighbor_[cellI]++;
            }
            else
            {
                //- Do nothing
            }
        }
    }

    //- Continue counting (for parallel run), update processor boundary
    forAll(filmPaperID_.boundaryField(),pI)
    {
        const fvPatchField<scalar>& paper_f = filmPaperID_.boundaryField()[pI];
    if(paper_f.type() == "processor")
    {
        //Info<<"patchType-processor: "<<paper_f.type()<<endl;
            const labelUList& faceCells = paper_f.patch().faceCells();
            vectorField faceB = paper_f.patch().Cf();
            forAll(paper_f, fi)
            {
                label cellO = faceCells[fi];
            scalar deltP = paper_f[fi] - filmPaperID_[cellO];
            if(deltP > 0)
            {
                filmPositiveDelNeighbor_[cellO]++;
                scalar dnbr = mag(Cv[cellO] - faceB[fi]);
                if(dnbr < filmMinPyrolysisDistance_[cellO])
                {
                    filmMinPyrolysisDistance_[cellO] = dnbr;
                    filmPyrolysisFrontCellPosition_[cellO] = faceB[fi];
                    filmPyrolysisFrontPaperID_[cellO] = filmPaperID_[cellO];
                }
            }
            else if(deltP < 0)
            {
                filmNegtiveDelNeighbor_[cellO]++;
            }
            else
            {
                //- Do nothing
            }
            }
    }
    }
    forAll(filmPaperID_, cellI)
    {
    if(filmNegtiveDelNeighbor_[cellI] > 0)
    {
        filmMinPyrolysisDistance_[cellI] = 1000;
        filmPyrolysisFrontCellPosition_[cellI] = vector(0,0,-1000.0);
        filmPyrolysisFrontPaperID_[cellI] = -filmPaperID_[cellI];
    }
    }
    
    // - update nearest pyrolysis front distance from neighbour cell
    for(int ipy=0;ipy<5;ipy++)
    {
        forAll(filmPyrolysisFrontCellPosition_,cellI)
        {
            //- For cells away from the pyrolysis front
        if(
        (filmPositiveDelNeighbor_[cellI] == 0) 
        && 
        (filmNegtiveDelNeighbor_[cellI] == 0)
          )
            {
                const labelList& nbrCL = regionMesh().cellCells()[cellI];
                forAll(nbrCL, nI)
                {
                    label nbrCellID = nbrCL[nI];
                    if(filmPaperID_[cellI] == filmPyrolysisFrontPaperID_[nbrCellID])
                    {
                        vector vnbr = Cv[cellI] - filmPyrolysisFrontCellPosition_[nbrCellID];
                        scalar dist = mag(vnbr);
                        if(dist < filmMinPyrolysisDistance_[cellI])
                        {
                            filmMinPyrolysisDistance_[cellI] = dist;
                            filmPyrolysisFrontCellPosition_[cellI] = 
                    filmPyrolysisFrontCellPosition_[nbrCellID];
                            filmPyrolysisFrontPaperID_[cellI] = filmPaperID_[cellI];
                        }
                    }
                }
            }
        }

        //- For parallel run, update processor boundary
        filmPyrolysisFrontCellPosition_.correctBoundaryConditions();
        filmPyrolysisFrontPaperID_.correctBoundaryConditions();

        forAll(filmPyrolysisFrontCellPosition_.boundaryField(),pI)
        {
            const fvPatchField<vector>& pyCpmin_f = 
            filmPyrolysisFrontCellPosition_.boundaryField()[pI];
    //if(pyCpmin_f.type() == "processor")
    //{
            const fvPatchField<scalar>& pyPage_f = 
            filmPyrolysisFrontPaperID_.boundaryField()[pI];
            const labelUList& faceCells = pyCpmin_f.patch().faceCells();
            vectorField faceB = pyCpmin_f.patch().Cf();
            forAll(pyCpmin_f, fi)
            {
                label cellO = faceCells[fi];
                //- For cells away from the pyrolysis front
            if(
                (filmPositiveDelNeighbor_[cellO] == 0) 
                && 
                (filmNegtiveDelNeighbor_[cellO] == 0)
              )
                {
                    vector vnbr = Cv[cellO] - pyCpmin_f[fi];
                    if(filmPaperID_[cellO] == pyPage_f[fi])
                    {
                        scalar dist = mag(vnbr);
                        if(dist < filmMinPyrolysisDistance_[cellO])
                        {
                            filmMinPyrolysisDistance_[cellO] = dist;
                            filmPyrolysisFrontCellPosition_[cellO] = pyCpmin_f[fi];
                            filmPyrolysisFrontPaperID_[cellO] = filmPaperID_[cellO];
                        }
                    }
                }
            }
        //}
    }
    }

    //- Correct bad cells 
    forAll(filmPaperID_, cellI)
    {
    if(fabs(filmPyrolysisFrontPaperID_[cellI]) != filmPaperID_[cellI])
        {
        Info<<"Correct-BAD-cells: "<<cellI<<tab<<filmPaperID_[cellI]<<tab
            <<filmPyrolysisFrontPaperID_[cellI]<<tab
        <<filmMinPyrolysisDistance_[cellI]<<tab
        <<filmPyrolysisFrontCellPosition_[cellI]<<endl;
            filmMinPyrolysisDistance_[cellI] = 1000;
            filmPyrolysisFrontCellPosition_[cellI] = vector(0,0,-1000.0);
        //filmPyrolysisFrontPaperID_[cellI] = filmPaperID_[cellI];
        //pyPage_[cellI] = paper_[cellI];
        }
    }

    //- Set thermally-thin zone based on minimum distance to pyFront
    forAll(filmMinPyrolysisDistance_,cellI)
    {
        scalar halfCellGridSize = 0.5*sqrt(cellVolume[cellI]*1000.0);
    scalar delPy = filmPeelingDistance_ - filmMinPyrolysisDistance_[cellI];
    if(delPy + halfCellGridSize < 0)
    {
        filmPaperPeelingZone_[cellI] = 0.0;
    }
    else if(delPy < halfCellGridSize)
    {
        filmPaperPeelingZone_[cellI] = (1.0 + delPy/halfCellGridSize)/2.0;
    }
    else
    {
        filmPaperPeelingZone_[cellI] = 1.0;
    }
    }

    //- Special treatment when pyrolysis fronts are stacked
    forAll(filmPaperID_, cellI)
    {
        if(
            (filmPositiveDelNeighbor_[cellI] > 0) 
            && 
            (filmNegtiveDelNeighbor_[cellI] > 0)
          )
        {
            filmPaperPeelingZone_[cellI] = 1.0;
        }

	    //- Turn on paper stack effect when paper sheet reach a threshold
	    if(
	        (filmPaperID_[cellI] > paperStackEffectID_) 
	        && 
	        (filmNegtiveDelNeighbor_[cellI] > 0)
	      )
	    {
	        filmPaperPeelingZone_[cellI] = 0.0;
	    }
    }

    // - Calculate maximum thermally thin region height
    const vectorField& cellCentre = regionMesh().cellCentres();
    filmPaperPyrolysisFront_ = dimensionedScalar("zero",dimLength,0.0);
    forAll(filmPaperPeelingZone_, cellI)
    {
        //scalar cellGridSize = sqrt(cellVolume[cellI]*1000.0);
        if(filmPaperPeelingZone_[cellI] > 0.5)
        {
            filmPaperPyrolysisFront_[cellI] = cellCentre[cellI].z();
        }
    }

    scalar maxPyHeight = max(filmPaperPyrolysisFront_).value();
    Info<<"ThinPaper-MaximumHeight: "<<time_.value()<<tab<<maxPyHeight<<endl;

    if(maxPyHeight < filmMinPyrolysisHeight_)
    {
        Info<<"Turnoff thermally thin model!"<<endl;
        forAll(filmPaperPeelingZone_, cellI)
        {
            filmPaperPeelingZone_[cellI] = 0;
        }
    }

    filmTpaper_.correctBoundaryConditions();
    filmTsurface_.correctBoundaryConditions();
    filmPaperPeelingZone_.correctBoundaryConditions();
    filmQnet_.correctBoundaryConditions();
    filmBlockFactor_.correctBoundaryConditions();
    filmMinPyrolysisDistance_.correctBoundaryConditions();
    filmMassReleaseRate_.correctBoundaryConditions();
    filmPositiveDelNeighbor_.correctBoundaryConditions();
    filmNegtiveDelNeighbor_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerRollPaper::thermoSingleLayerRollPaper
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    thermoSingleLayer(modelType, mesh, g, regionType, readFields),
    filmPaperID_
    (
        IOobject
        (
            "paperID",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmTpaper_
    (
        IOobject
        (
            "paperTemperature",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmTsurface_
    (
        IOobject
        (
            "paperSurfaceTemperature",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 300.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    filmPaperPeelingZone_
    (
        IOobject
        (
            "peelingZone",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmQnet_
    (
        IOobject
        (
            "Qnet",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmBlockFactor_
    (
        IOobject
        (
            "blockFactor",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmMassReleaseRate_
    (
        IOobject
        (
            "MLR",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmPaperPyrolysisFront_
    (
        IOobject
        (
            "pyrolysisFront",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    //- new model variable to calculate distance
    filmPositiveDelNeighbor_
    (
        IOobject
        (
            "positiveDelNeighbor",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmNegtiveDelNeighbor_
    (
        IOobject
        (
            "negtiveDelNeighbor",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmMinPyrolysisHeight_(coeffs_.lookupOrDefault<scalar>("minPyHeight", 0.0)),
    filmPeelingDistance_(coeffs_.lookupOrDefault<scalar>("peelingDistance", 0.05)),
    paperStackEffectID_(coeffs_.lookupOrDefault<scalar>("paperStack", GREAT)),

    filmMinPyrolysisDistance_
    (
        IOobject
        (
            "minPyDistance",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimLength, 1000.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmPyrolysisFrontCellPosition_
    (
        IOobject
        (
            "minPyPosition",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedVector("zero", dimLength, vector(0,0,-1000)),
        zeroGradientFvPatchScalarField::typeName
    ),

    filmPyrolysisFrontPaperID_
    (
        IOobject
        (
            "pyFrontPaperID",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )



{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerRollPaper::~thermoSingleLayerRollPaper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
