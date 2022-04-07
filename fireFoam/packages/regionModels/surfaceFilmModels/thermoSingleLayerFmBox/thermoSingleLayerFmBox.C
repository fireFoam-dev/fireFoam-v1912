/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "thermoSingleLayerFmBox.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedWallPolyPatch.H"
#include "mathematicalConstants.H" 
#include "cyclicFvPatch.H"


// Sub-models
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"
#include "stdio.h"
#include "assert.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayerFmBox, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerFmBox, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayerFmBox::read()
{
    if(thermoSingleLayerPw::read())
    {
        trappedMassModel_ = coeffs_.lookupOrDefault<Switch>("trappedMassModel",false);
        if(trappedMassModel_){
            const dictionary& subdict=coeffs_.subDict("trappedMassModelCoeffs");
            subdict.lookup("fraction") >> trapFraction_; 
            subdict.lookup("profileFactor") >> trapProfileFactor_; 
            subdict.lookup("trapMax") >> trapMax_; 
            subdict.lookup("boxes") >> nBoxes_; 
            subdict.lookup("tiers") >> nTiers_; 
        }
        return true;
    }
    else
    {
        return false;
    }
}



void thermoSingleLayerFmBox::assignBoxTierSide(){ 
    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();
    forAll(bm,patchi){
        const labelUList& bCells = bm[patchi].faceCells();
//        const labelUList& nbrFaceCells =
//                bm[patchi].neighbPatch().faceCells(); //gets cyclic patch neighbor face cells

//        const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
//        (
//            bm[patchi]
//        ).neighbFvPatch();
//        const labelList& nbrFaceCells = nbrPatch.patch().faceCells();

        const label sides=6;
        for(int boxIndex=1;boxIndex<=nBoxes_;boxIndex++){
            char buffer[80];
            sprintf(buffer,"%02d",boxIndex);
            string cubeName="cube";
            cubeName.append(buffer);


            for(int sideIndex=1;sideIndex<=sides;sideIndex++){

//            	match "top" boundaries with string "cube01_side1_top" and assign box, side, and tier indices
            	sprintf(buffer,"%01d",sideIndex);
                string sideName="side";
                sideName.append(buffer);
                string name;
                name.append(cubeName);
                name.append("_");
                name.append(sideName);
                name.append("_top");
                //Info << name.c_str()<<endl;
                if(bm[patchi].name()==name){ //top tier
              //  	Info << "setting box, side, and tier for " << name << endl;
                	forAll(bCells,facei){
                        const label index=bCells[facei];
                        box_[index]=boxIndex;
                        side_[index]=sideIndex;
                        tier_[index]=boxIndex;
                    }
                }


  /*            find box-top boundaries (sideIndex==5) and define edge faces as topEdge
                cube01_side5_to_cube01_side1 means the cyclic boundary between side5 and side1 (from the perspective of side5)
                cube01_side1_to_cube01_side5 means the cyclic boundary between side1 and side5 (from the perspective of side1)*/
                {
                	string edgeName;
                	edgeName.append(cubeName);
                	edgeName.append("_");
                	label sideIndex2=5; //top surface
                	sprintf(buffer,"%01d",sideIndex2);
                	string sideName2="side";
                	sideName2.append(buffer);
                	edgeName.append(sideName2);
                	edgeName.append("_to_");
                	edgeName.append(cubeName);
                	edgeName.append("_");
                	string name;
                	name.append(edgeName);
                	name.append(sideName);
                	if(bm[patchi].name()==name)
                	{
                		forAll(bCells,facei){
                			const label index=bCells[facei];
                			topEdge_[index]=1;
                		}
                	}
                }
                //find box-bottom boundaries (sideIndex==6) and define edge faces as bottomEdge
                {
                	string edgeName;
                	edgeName.append(cubeName);
                	edgeName.append("_");
                	label sideIndex2=6; //=bottom surface
                	sprintf(buffer,"%01d",sideIndex2);
                	string sideName2="side";
                	sideName2.append(buffer);
                	edgeName.append(sideName2);
                	edgeName.append("_to_");
                	edgeName.append(cubeName);
                	edgeName.append("_");
                	string name;
                	name.append(edgeName);
                	name.append(sideName);
                	if(bm[patchi].name()==name)
                	{
                		forAll(bCells,facei){
                			const label index=bCells[facei];
                			bottomEdge_[index]=1;
                		}
                	}
                }
            }
        }
    }
//    assert(1==0);
}

void thermoSingleLayerFmBox::integrateMass(){ 
        const label sides=6;

        static DynamicList<scalar> massBoxes;

        massBoxes.resize(nBoxes_);
        forAll(massBoxes,i){
            massBoxes[i]=0;
        }
        forAll(delta_,index){
            massBoxes[static_cast<int>(box_[index]-1)]+=(delta_[index])*rho_[index]*magSf()[index];

        }
        //TODO:  find better 'one-liner' to reduce a fixedlist
        forAll(massBoxes,i){
            reduce(massBoxes[i], sumOp<scalar>());
        }
        if (debug){
            char buffer[256];
            sprintf(buffer,"time = %10.5f wholeBox ",time_.time().value());
            Info << buffer;
            for(int boxIndex=1;boxIndex<=nBoxes_;boxIndex++){
                sprintf(buffer,"\t %2d = %10.5g \t",boxIndex,massBoxes[boxIndex-1]);
                Info << buffer;            
            }
            Info << endl;
        }

        /*sides*/
        static DynamicList<FixedList<scalar, sides> > massSides;
        massSides.resize(nBoxes_);
        forAll(massSides,i){
            massSides[i]=0.0;
        }

        forAll(delta_,index){
            massSides[static_cast<label>(box_[index]-1)][static_cast<label>(side_[index]-1)]+=(delta_[index])*rho_[index]*magSf()[index];
        }
        //TODO:  find better 'one-liner' to reduce a fixedlist
        forAll(massSides,i){
            forAll(massSides[i],j){
                reduce(massSides[i][j], sumOp<scalar>());
            }
        }
        if (debug){
            char buffer[256];
            for(int boxIndex=1;boxIndex<=nBoxes_;boxIndex++){
                sprintf(buffer,"time = %10.5f box%1dsides ",time_.time().value(),boxIndex);
                Info << buffer;
                for(int sideIndex=1;sideIndex<=sides;sideIndex++){
                    sprintf(buffer,"%2d,%1d = %10.5g ",boxIndex,sideIndex,massSides[boxIndex-1][sideIndex-1]);
                    Info << buffer;            
                }
                Info << endl;
            }
        }
}

void thermoSingleLayerFmBox::initialise()
{
    if (debug)
    {
        Pout<< "thermoSingleLayerFmBox::initialise()" << endl;
    }

    read();
    if(trappedMassModel_){
        assignBoxTierSide();
    }
}

void thermoSingleLayerFmBox::updateSubmodels()
{
    thermoSingleLayerPw::updateSubmodels();

    if(trappedMassModel_){
        updateMassLost();
        integrateMass();
        rhoSpLost_=-massLost_/magSf()/time_.deltaT();
        rhoSp_+=rhoSpLost_;
    }
}

void thermoSingleLayerFmBox::updateMassLost()
{
    const dimensionedScalar deltaMin("deltaMin",dimLength,1e-10); 
    const scalarField availableMass = max(dimensionedScalar("zero", dimMass, 0.0),(delta_ - deltaMin)*rho_*magSf());

    trapSum.resize(nTiers_);
    forAll(trapSum,i){
        trapSum[i]=0;
    }
    forAll(tier_,index){
        if(side_[index]==5&&topEdge_[index]==1){
            trapSum[static_cast<int>(tier_[index]-1)]+=trappedMass_[index];
        }
    }
    //TODO:  find better 'one-liner' to reduce a fixedlist
    forAll(trapSum,i){
        reduce(trapSum[i], sumOp<scalar>());
    }
    
    forAll(tier_,index){
        label tierIndex=static_cast<int>(tier_[index]);
        if(tierIndex==1){
            if(side_[index]==5){
                /*TODO: next line is the same as mf=trapSum[tierIndex]-1]/(trapMax.value()+SMALL);*/
                const scalar mf=1.0-(trapMax_.value()-trapSum[tierIndex-1])/(trapMax_.value()+SMALL);
                scalar trapScale=1.0-(exp(trapProfileFactor_*mf)-1.0)/(exp(trapProfileFactor_+ROOTVSMALL)-1.0);
                scalar trap=0.0;
                if(topEdge_[index]==1){
                    if(trapSum[tierIndex-1]<trapMax_.value()){
                        trap=trapScale*trapFraction_*availableMass[index];
                    }
                }
                trappedMass_[index]+=trap;
                massLost_[index]=trap;
            }
        }
        else if(tierIndex==2){
            if(side_[index]==5){
                const scalar mf=1.0-(trapMax_.value()-trapSum[tierIndex-1])/(trapMax_.value()+SMALL);
                scalar trapScale=1.0-(exp(trapProfileFactor_*mf)-1.0)/(exp(trapProfileFactor_)-1.0);
                scalar trap=0.0;
                if(topEdge_[index]==1){
                    if(trapSum[tierIndex-1]<trapMax_.value()){
                        trap=trapScale*trapFraction_*availableMass[index];
                    }
                }
                trappedMass_[index]+=trap;
                massLost_[index]=trap;
            }
        }
        else if(tierIndex==3){
            if(side_[index]==5){
                const scalar mf=1.0-(trapMax_.value()-trapSum[tierIndex-1])/(trapMax_.value()+SMALL);
                scalar trapScale=1.0-(exp(trapProfileFactor_*mf)-1.0)/(exp(trapProfileFactor_)-1.0);
                scalar trap=0.0;
                if(topEdge_[index]==1){
                    if(trapSum[tierIndex-1]<trapMax_.value()){
                        trap=trapScale*trapFraction_*availableMass[index];
                    }
                }
                trappedMass_[index]+=trap;
                massLost_[index]=trap;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerFmBox::thermoSingleLayerFmBox
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    thermoSingleLayerPw(modelType, mesh, g, regionType, false),
    trappedMassModel_(false),
    nBoxes_(0),
    nTiers_(0),
    trapFraction_(0.0), 
    trapProfileFactor_(0.0), 
    trapMax_("trapMax",dimMass,0.0), 
    trappedMass_
    (
        IOobject
        (
            "trappedMass",
            time_.timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    massLost_
    (
        IOobject
        (
            "massLost",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    rhoSpLost_
    (
        IOobject
        (
            "massLostSourceTerm",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0)
    ),
    box_
    (
        IOobject
        (
            "box",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
    ),
    tier_
    (
        IOobject
        (
            "tier",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
    ),
    side_
    (
        IOobject
        (
            "side",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
    ),
    topEdge_
    (
        IOobject
        (
            "topEdge",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
    ),
    bottomEdge_
    (
        IOobject
        (
            "bottomEdge",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
         0.0,
         zeroGradientFvPatchScalarField::typeName
    )
    //outFilmExtent_("postProcessing/filmExtent/filmExtent_"+time_.timeName())
{
    initialise();

    //outFilmExtent_ << "# Time[s]    " << "filmExtentMin(x,y,z)  " << "filmExtentMax(x,y,z)  "  << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerFmBox::~thermoSingleLayerFmBox()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayerFmBox::info()
{
    thermoSingleLayerPw::info();

    static char buffer[256];
    forAll(trapSum,i){
        sprintf(buffer,"trapSum[%d]         = %-10.5g\n",i+1,trapSum[i]);
        Info << indent<<buffer;
    }
    //cout.precision(10);
    //Info << indent;
    //cout.width(5);
    //Info << endl;
    //Info << time_.time().value() << " ";
    //forAll(trapSum,i){
    //    Info << "trapSum["<<i<<"] = "<<trapSum[i]<< " ";
    //}
    //Info << endl;
    
    // film extent calculation
    //static label count=0;
    //if(count%10==0)
    //{
    //     
    //    dimensionedVector filmExtentMax("zero",dimLength,Foam::vector(0,0,0));
    //    dimensionedVector filmExtentMin("zero",dimLength,Foam::vector(0,0,0));

    //    const pointField& cellCentres = regionMesh().cellCentres();

    //    forAll (omega_.internalField(),cellI)
    //    {
    //        if (omega_[cellI] >= 1)
    //        {
    //            if (cellCentres[cellI][0] > filmExtentMax[0].value()) 
    //            {
    //                DEBUG("filmExtentMax set");
    //                //filmExtentMax[0].value() = dimensionedScalar("val",dimLength,cellCentres[cellI][0]);
    //                filmExtentMax[0].value() = cellCentres[cellI][0];
    //            }
    //            DEBUG(cellCentres[cellI][0]);
    //            DEBUG(filmExtentMax[0].value());
    //            filmExtentMax[1] = dimensionedScalar("val",dimLength,min( cellCentres[cellI][1], filmExtentMax[1].value()));
    //            filmExtentMax[2] = dimensionedScalar("val",dimLength,max( cellCentres[cellI][2], filmExtentMax[2].value()));
    //            filmExtentMin[0] = dimensionedScalar("val",dimLength,min( cellCentres[cellI][0], filmExtentMin[0].value()));
    //            filmExtentMin[1] = dimensionedScalar("val",dimLength,min( cellCentres[cellI][1], filmExtentMin[1].value()));
    //            filmExtentMin[2] = dimensionedScalar("val",dimLength,min( cellCentres[cellI][2], filmExtentMin[2].value()));
    //        }
    //    }

    //    outFilmExtent_ <<  time_.value() << "  "
    //        << returnReduce(filmExtentMax[0].value(),maxOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMax[1].value(),maxOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMax[2].value(),maxOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMin[0].value(),minOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMin[1].value(),minOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMin[2].value(),minOp<scalar>()) << "  "
    //        << endl;
    //    Info <<  time_.value() << "  "
    //        << returnReduce(filmExtentMax[0].value(),maxOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMax[1].value(),maxOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMax[2].value(),maxOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMin[0].value(),minOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMin[1].value(),minOp<scalar>()) << "  "
    //        << returnReduce(filmExtentMin[2].value(),minOp<scalar>()) << "  "
    //        << endl;
    //}
    //count++;


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam


// ************************************************************************* //
