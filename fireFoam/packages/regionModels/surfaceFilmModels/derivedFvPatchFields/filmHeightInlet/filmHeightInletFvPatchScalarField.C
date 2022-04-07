/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "filmHeightInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Random.H" //kvm
#include "constants.H" //kvm

#include <iostream>
#include <fstream>
using namespace std;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmHeightInletFvPatchScalarField::
filmHeightInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    deltafName_("deltaf"),
    mean_(0.0002),
    amplitude_(0.00002),
    spatialFrequency_(1.0),
    temporalFrequency1_(1.0),
    temporalFrequency2_(1.0)
{}


Foam::filmHeightInletFvPatchScalarField::
filmHeightInletFvPatchScalarField
(
    const filmHeightInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    deltafName_(ptf.deltafName_),
    mean_(ptf.mean_),
    amplitude_(ptf.amplitude_),
    spatialFrequency_(ptf.spatialFrequency_),
    temporalFrequency1_(ptf.temporalFrequency1_),
    temporalFrequency2_(ptf.temporalFrequency2_)
{}


Foam::filmHeightInletFvPatchScalarField::
filmHeightInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    deltafName_(dict.lookupOrDefault<word>("deltaf", "deltaf")),
    mean_(dict.lookupOrDefault<scalar>("mean", 0.0002)),
    amplitude_(dict.lookupOrDefault<scalar>("amplitude", 0.00002)),
    spatialFrequency_(dict.lookupOrDefault<scalar>("spatialFrequency", 1)),
    temporalFrequency1_(dict.lookupOrDefault<scalar>("temporalFrequency1", 1)),
    temporalFrequency2_(dict.lookupOrDefault<scalar>("temporalFrequency2", 1))
{
    /*Info << "reading mean\n";*/
    /*Info << "mean " << mean_<<endl;*/
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::filmHeightInletFvPatchScalarField::
filmHeightInletFvPatchScalarField
(
    const filmHeightInletFvPatchScalarField& fhivpvf
)
:
    fixedValueFvPatchScalarField(fhivpvf),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_),
    mean_(fhivpvf.mean_),
    amplitude_(fhivpvf.amplitude_),
    spatialFrequency_(fhivpvf.spatialFrequency_),
    temporalFrequency1_(fhivpvf.temporalFrequency1_),
    temporalFrequency2_(fhivpvf.temporalFrequency2_)
{}


Foam::filmHeightInletFvPatchScalarField::
filmHeightInletFvPatchScalarField
(
    const filmHeightInletFvPatchScalarField& fhivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fhivpvf, iF),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_),
    mean_(fhivpvf.mean_),
    amplitude_(fhivpvf.amplitude_),
    spatialFrequency_(fhivpvf.spatialFrequency_),
    temporalFrequency1_(fhivpvf.temporalFrequency1_),
    temporalFrequency2_(fhivpvf.temporalFrequency2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmHeightInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    /*const fvsPatchField<scalar>& phip =*/
        /*patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);*/

    /*const fvPatchField<scalar>& rhop =*/
        /*patch().lookupPatchField<volScalarField, scalar>(rhoName_);*/

    const fvPatchField<scalar>& deltafp =
        patch().lookupPatchField<volScalarField, scalar>(deltafName_);

    scalarField dfa(deltafp.size()); 
    /*scalarList dfa(deltafp.size()); */
    /*fvPatchField<scalar> dfa; */
    /*dfa=deltafp;*/

    /*vectorField n = patch().nf();*/
    /*const scalarField& magSf = patch().magSf();*/
    /*const volVectorField& cellCentres = this->patch().patch().faceCentres();*/
    const pointField& centers = this->patch().patch().faceCentres();

    /*operator==(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + ROOTVSMALL));*/
    scalar time=this->db().time().timeOutputValue();
    /*scalar freq=1;*/
    /*scalar amp=0.00001;*/
    /*scalar df=mean_+amp*sin(time*freq*2.0*constant::mathematical::pi);*/
    const scalar twoPi=2.0*constant::mathematical::pi;
    // const scalar pi=constant::mathematical::pi;
    scalar tempOsc=cos(time*temporalFrequency1_*twoPi);
    scalar tempOsc2=cos(time*temporalFrequency2_*twoPi);
    /*Info << "time " << time << " df " <<df<<endl;*/
    /*static label seed=0;*/
    /*seed=seed+1;*/
    /*static Random ranGen_(seed);//kvm*/

    static label direction=0;
    static scalar length=0;
    static bool first=true;
    if(first){
        getDirectionLength(direction,length);
        first=false;
    }

    /*scalar freqSpatial=1.0;*/
    /*scalar spatOsc=sin(time*freqSpatial*twoPi);*/
    //scalar meanL=0.00018; 
    /*scalar meanL=0.0006; */
    /*scalar amplitudeL=0.0000200; // m/s*/
    /*scalar peaks=4;*/
    /*scalar period=length/2/pi;*/
    /*scalar frequency=peaks/period;*/
    forAll(dfa,i){
        /*scalar rand=0.0;*/
        /*ranGen_.randomise(rand); */
        /*rand-=0.5; */
        /*Info << "rand " << rand<<endl;*/

        dfa[i]=mean_+tempOsc*amplitude_*mean_*sin(spatialFrequency_*twoPi*centers[i][direction]+tempOsc2*twoPi);
        /*dfa[i]=mean_+tempOsc*amplitude_*sin(spatialFrequency_*twoPi*centers[i][direction]+temporalFrequency2_*twoPi*time);*/
        /*dfa[i]=mean_+amplitude_*sin(spatialFrequency_*twoPi*centers[i][direction]+temporalFrequency2_*twoPi*time);*/
        /*dfa[i]=df+0.0001*rand;*/
    }
    /*operator==(mean_ + 0.0001*rand);*/
    /*Info << "dfa"<<endl;*/

    if(0){
        ofstream myfile;
        char buffer[256];
        sprintf(buffer,"dfa_%010.5f",time);
        myfile.open (buffer);
        forAll(dfa,i){
            myfile << dfa[i] <<" \n";
        }
        myfile.close();
        /*Info << dfa;*/
    }
    operator==(dfa);

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::filmHeightInletFvPatchScalarField::getDirectionLength(label &direction,scalar& length)
{
    scalar minX=1e10;
    scalar maxX=-1e10;
    scalar minY=1e10;
    scalar maxY=-1e10;
    scalar minZ=1e10;
    scalar maxZ=-1e10;
    const pointField& centers = this->patch().patch().faceCentres();
    forAll(centers,i){
        //find xLength
        if(centers[i][0]<minX){
            minX=centers[i][0];
        }
        if(centers[i][0]>maxX){
            maxX=centers[i][0];
        }
        //find yLength
        if(centers[i][1]<minY){
            minY=centers[i][1];
        }
        if(centers[i][1]>maxY){
            maxY=centers[i][1];
        }
        //find zLength
        if(centers[i][2]<minZ){
            minZ=centers[i][2];
        }
        if(centers[i][2]>maxZ){
            maxZ=centers[i][2];
        }
    }
    scalar lengthX=maxX-minX;
    scalar lengthY=maxY-minY;
    scalar lengthZ=maxZ-minZ;
    if(lengthX>lengthY&&lengthX>lengthZ){
        direction=0;
        length=lengthX;
    }
    if(lengthY>lengthX&&lengthY>lengthZ){
        direction=1;
        length=lengthY;
    }
    if(lengthZ>lengthX&&lengthZ>lengthY){
        direction=2;
        length=lengthZ;
    }
    return;
}

void Foam::filmHeightInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("deltaf", "deltaf", deltafName_);
    /*mean_.writeEntry("mean",os);*/
    os.writeKeyword("mean") << mean_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("spatialFrequency") << spatialFrequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("temporalFrequency1") << temporalFrequency1_ << token::END_STATEMENT << nl;
    os.writeKeyword("temporalFrequency2") << temporalFrequency2_ << token::END_STATEMENT << nl;
    /*writeEntryIfDifferent<scalar>(os, "mean", 0.0002, deltafName_);*/
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::filmHeightInletFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& pvf
)
{
    fvPatchField<scalar>::operator=(pvf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        filmHeightInletFvPatchScalarField
    );
}


// ************************************************************************* //
