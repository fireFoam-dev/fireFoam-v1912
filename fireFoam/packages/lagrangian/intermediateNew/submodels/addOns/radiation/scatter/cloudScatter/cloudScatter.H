/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::radiation::cloudScatter

Description
    Cloud radiation scatter model

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_cloudScatter_H
#define radiation_cloudScatter_H

#include "scatterModel.H"
#include "fvDOM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
// Forward declaration of classes
    class fvDOM; // luwi
/*---------------------------------------------------------------------------*\
                         Class cloudScatter Declaration
\*---------------------------------------------------------------------------*/

class cloudScatter
:
    public scatterModel
{

public:

//<fmglobal>
    enum phaseFunctionType
    {
        preCorrected,
        henyeyGreenstein,
        hgDirect,
	hgfixedAsymmetry
    };

    // Private data
private:
        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Cloud names
        wordList cloudNames_;

       //- Dictionary keywords for specifying different phase function options
       static Enum<phaseFunctionType> phaseFunctionTypeNames_; 

       //- Enumerator that identifies phase function option
       const phaseFunctionType phaseFoption_; 

       //- Fixed value of asymmetry factor for: phaseFoption_ = hgfixedAsymmetry
       const scalar asymmetryFactor_;

       //- Minimum solid angle resolution for averaging phase function
       const label nminTheta_;
       const label nminPhi_;
    
       //- Number of fvDOM rays
       label nRays_; 

       //- Number of radiation spectral bands
       label nBands_; 

       //- Number of unique scattering phase function values
       label nphaseFUniq_; 

       //- Number of unique values of phase function normalization factor: phaseF[][j]*omega[]
       label ncscatUniq_; 

       //- Row indices of the unique phase function values
       labelList iphaseFUniqRow_; 

       //- Column indices of the unique phase function values
       labelList iphaseFUniqCol_; 

       //- Row indices of the unique normalization factors
       List<label> icscatUniq_; 

       //- row indexes of phase function normalization factor in the reduced list "cscatUniq_"
       labelList cscatRow_; 

       //- Map of phase function unique values into the full matrix
       List<labelList> phaseFindx_; 

       //- phaseFUniq_[i][j] contains volScalarField of the jth unique phase function for lambda i
       PtrList<PtrList<volScalarField>> phaseFUniq_; 

       //- phaseFUniq_ at a reference value of aymmetry factor (used for hgfixedAsymmetry)
       scalarList refPhaseFUniq_;
    
       //- cosines of angles between direction vectors
       scalarList cosTHETA_;

       //- Unique normalization factor values
       PtrList<PtrList<volScalarField>> cscatUniq_; 

       //- cscatUniq_ at a reference value of aymmetry factor (used for hgfixedAsymmetry)
       scalarList refCscatUniq_;    

       //- Switch to indicate whether arrays phaseFUniq_ and cscatUniq_ need to be initialized
       Switch doInit_; 
//</fmglobal>
    
public:

    //- Runtime type information
    TypeName("cloudScatter");


    // Constructors

        //- Construct from components
        cloudScatter
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cloudScatter();


    // Member Operators

        // Access

            //- Return scatter coefficient
//<fmglobal>    
            virtual tmp<volScalarField> sigmaEff(const label bandI, const label iRay, const scalar omega) const; 

            //- Return phase-function for computing in-scattering 
            virtual tmp<volScalarField> pFunc(const label bandI, const label sour, const label dest) const; 

        // Edit
            //- Setup phaseFUniq_ and cscatUniq_; and set flag to indicate this is done
            void initPhaseFuncs(const fvDOM& dom); 

            //- Setup arrays for in-scattering calculations
            void setupMaps(const fvDOM& dom); 

            //- Pre-compute (fixed) scattering phase functions averaged over solid angles
            void angleAveragedPhaseF(const fvDOM& dom);     

            //- Update phaseFUniq_ and cscatUniq_
            virtual void updatePhaseFuncs(const fvDOM& dom); 

            //- Flag to indicate cloud model is being used with PRECORRECTED
            virtual Switch cloudCorrectedSigma() const; 
//</fmglobal>
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
