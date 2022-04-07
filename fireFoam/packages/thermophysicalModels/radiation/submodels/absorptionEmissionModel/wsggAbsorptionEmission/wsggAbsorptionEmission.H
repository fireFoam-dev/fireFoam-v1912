/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::radiation::wsggAbsorptionEmission

Description

    Johansson:
          Account for variations in the H2O to CO2 molar ratio when modelling gaseous
      radiative heat transfer with the weighted-sum-of-grey-gases model
      Robert Johansson (corr author), Bo Leckner, Klas Andersson, Filip Johnsson
      Combustion and Flame 158(2011) 893-901
    Franca:
          WSGG correlations based on HITEMP2010 for computation of thermal
      radiation in non-isothermal, non-homogeneous H2O/CO2 mixtures
      Leonardo J. Dorigon, Gustavo Duciak, Rogerio Brittes, Fabiano Cassol,
      Marcelo Galarca, and Francis H.R. Franca
      International Journal of Heat and Mass Transfer 64(2013) 863-873
    Smith:
      Evaluation of Coefficients for the Weighted Sum of Grey Gases Model
      T.F. Smith, Z.F. Shen, J.N. Friedman
      Journal of Heat Transfer, Vol. 104 (Nov. 1982), pgs 602-608


    wsgg dictionary (must be in file: constant/wsggProperties):

    \verbatim
    johansson
    {
       CandF2011
       {
          Tref 1200;
      grayGas1
      {
         Tcommon         298.;
         Tlow            298.;
         Thigh           2500.;
         invTemp         false;
         loTcoeffs ();
         hiTcoeffs (0.055 0.012 0.358 0.0731 -0.0466 -0.165 -0.0554 0.093 0.0598 0.0028 -0.0256);
      }
      greyGas2
      {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (0.88 -0.021 0.392 -0.212 0.0191 -0.291 0.644 -0.209 0.0784 -0.197 0.0662);
      }
      greyGas3
      {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (10 -1.6 0.142 -0.0831 0.0148 0.348 -0.294 0.0662 -0.122 0.118 -0.0295);
      }
      grayGas4
      {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (135 -35 0.0798 -0.037 0.0023 0.0866 -0.106 0.0305 -0.0127 0.0169 -0.0051);
      }
      }
    }

    franca
    {
      mixMoleRatio2
      {
    grayGas1
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.921e-1 5.617e-2 7.844e-4 -8.563e-7 4.246e-10 -7.440e-14);
    }
    greyGas2
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.719 1.426e-1 1.795e-4 -1.077e-8 -6.971e-11 1.774e-14);
    }
    greyGas3
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.137e1 1.362e-1 2.574e-4 -3.711e-7 1.575e-10 -2.267e-14);
    }
    grayGas4
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.110e2 1.222e-1 -2.327e-5 -7.492e-8 4.275e-11 -6.608e-15);
    }
      }

      mixMoleRatio1
      {
    grayGas1
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.873e-1 7.197e-2 8.724e-4 -9.690e-7 4.651e-10 -7.917e-14);
    }
    greyGas2
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.723 1.107e-1 3.397e-4 -2.467e-7 4.647e-11 -1.039e-15);
    }
    greyGas3
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.248e1 2.091e-1 -6.423e-5 -3.200e-8 1.718e-11 -2.105e-15);
    }
    grayGas4
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2500.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (1.449e2 7.092e-2 6.586e-5 -1.278e-7 5.577e-11 -7.709e-15);
    }
      }

    }

    smith
    {
      mixMoleRatio2
      {
    grayGas1
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2400.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (0.4201 6.508e-1 -5.551e-4 3.029e-7 -5.353e-11);
    }
    greyGas2
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2400.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (6.516 -0.2504e-1 6.112e-4 -3.882e-7 6.528e-11);
    }
    greyGas3
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2400.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (131.9 2.718e-1 -3.118e-4 1.221e-7 -1.612e-11);
    }
      }

      mixMoleRatio1
      {
    grayGas1
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2400.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (0.4303 5.150e-1 -2.303e-4 0.9779e-7 -1.494e-11);
    }
    greyGas2
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2400.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (7.055 0.7749e-1 3.399e-4 -2.297e-7 3.770e-11);
    }
    greyGas3
    {
        Tcommon         298.;
        Tlow            298.;
        Thigh           2400.;
        invTemp         false;
        loTcoeffs ();
        hiTcoeffs (178.1 1.907e-1 -1.824e-4 0.5608e-7 -0.5122e-11);
    }
      }

    }

    \endverbatim


SourceFiles
    wsggAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef wsggAbsorptionEmission_H
#define wsggAbsorptionEmission_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "fluidThermo.H"

#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "reactingMixture.H"
#include "blackBodyEmission.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of bands
        static const int maxBands_ = 10;

        //-  Absorption coefficients
        FixedList<absorptionCoeffs, maxBands_> coeffs_;

       //- WSGG model selection
       enum correlationType
       {
           absCoeffJohansson,
           absCoeffFranca,
           absCoeffSmith
       };
    
private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Choice of available WSGG correlation types
        const word wsggModelName_;

        //- Choice of available WSGG correlations of the selected type
        word correlationID_;
    
        //- Enumerator that identifies phase function option
        correlationType wsggModel_;

        //- Name of file containing the wsgg dictionary
        const word wsggFileName_;

        //- wsgg dictionary
        const dictionary wsggDict_;
    
        //- Thermo package
        const reactingMixture<gasHThermoPhysics>& gasThermo_;

        //- pressure
        const volScalarField& p_;
        //- temperature
        const volScalarField& T_;

        //- Species mass fractions
        const PtrList<volScalarField>& Y_;
    
        //- Species index of CO2 in mixture
        const label indexCO2_;

        //- Species index of H2O in mixture
        const label indexH2O_;

        //- Molecular weight of CO2
        const scalar molWtCO2_;

        //- Molecular weight of H2O
        const scalar molWtH2O_;
    
        //- Bands
        label nBands_;

        //- Access to soot coefficient from case file
        const scalar Csoot_;

        //- Flag to indicate whether soot model is mixtureFractionSoot
        const Switch mixtureFractionSoot_;

        //- Get thermo refs, indexes and values
        const psiReactionThermo& thermoRef(const fvMesh& mesh) const;

        //- Get wsgg dictionary
        const dictionary wsggDict(const fvMesh& mesh);
    
        //- Using WSGG in gray or non-gray simulation?
        const Switch isGray_;

        //- Read user specification of spectralMethod (gray or banded)
        const Switch getGrayFlag(const dictionary& dict) const;

        //- Characteristic path length / [m] from dict
        const dimensionedScalar pathLength_;

        //- Reference temperature for normalization in Johansson fit of WSGG weights
        const dimensionedScalar Tref_;

       //- Pressure conversion factors
       const scalar paToBar_ = 1.e-5;
       const scalar paToAtm_ = 9.869231e-6;

public:

    //- Runtime type information
    TypeName("wsggAbsorptionEmission");


    // Constructors

        //- Construct from components
        wsggAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;

        // weighting coefficient for wsggm
                tmp<volScalarField> ggWeight(const label bandI = 0) const;

            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


        inline bool isGrey() const
        {
            return isGray_;
        }

        //- Number of bands
        inline label nBands() const
        {
            return ( isGray_? 1 : nBands_ );
        }

        void correct
        (
            volScalarField& a_,
            PtrList<volScalarField>& aLambda

        ) const;

        virtual void correctEnFrac
        (
            PtrList<volScalarField>& enFrac,
            blackBodyEmission& blackBody
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
