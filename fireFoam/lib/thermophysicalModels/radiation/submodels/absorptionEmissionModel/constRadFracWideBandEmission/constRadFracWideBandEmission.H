/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::radiation::constRadFracWideBandEmission

Description
    constRadFracWideBandEmission radiation absorption and emission coefficients
    for continuous phase

    There is no check of continuity of the bands. They should not ovelap or
    have gaps.

   band dictionary:
   \verbatim
       band0
       {
           bandLimits (1.0e-6 2.63e-6);
       }
   \endverbatim


SourceFiles
    constRadFracWideBandEmission.C

\*---------------------------------------------------------------------------*/

#ifndef constRadFracWideBandEmission_H
#define constRadFracWideBandEmission_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "basicThermo.H"

//#include "OpenFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class constRadFracWideBandEmission Declaration
\*---------------------------------------------------------------------------*/

class constRadFracWideBandEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        // Maximum number of species considered for absorptivity
        static const int nSpecies_ = 5;

        // Absorption Coefficients
        absorptionCoeffs coeffs_[nSpecies_];


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Bands
        List<Vector2D<scalar>> iBands_;

        //- Proportion of the energy emitted for various bands
        scalarList enFracs_;

        //- Emission constant coefficient
        const scalar EhrrCoeff_;

        //- Switch for scaling radiant fraction
        Switch radScaling;

        //- Radiant fraction for fuel from patch 1 and 2
        scalar Ehrr1_;
        scalar Ehrr2_;

        //- patch names for patches 1 and 2
        wordList patchName1_;
        wordList patchName2_;

        //- min and max specified fuel radiant fractions
        const scalar EhrrMin_;
        const scalar EhrrMax_;

        //- Radiation temperature
        scalar radTemp_;

        //- Bands
        label nBands_;
        
        // Bool for energy fractions for various bands..
        bool fracsSet_;

        // function to set energy fractions for various bands
        void setEnergyFracs();

public:

    //- Runtime type information
    TypeName("constRadFracWideBandEmission");


    // Constructors

        //- Construct from components
        constRadFracWideBandEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constRadFracWideBandEmission();


    // Member Operators

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI) const;


    // Member Functions

        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }

        //- Lower and upper limit of band i
        inline const Vector2D<scalar>& bands(const label i) const
        {
            return iBands_[i];
        }

        //- Update routine
        void correct
        (
            volScalarField& a_,
            PtrList<volScalarField>& aLambda
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
