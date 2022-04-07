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
    Foam::radiation::boxFixedBandsAbsorptionEmission

Description

    Box model for gas radiation, using fixed spectral bands. This is a gray-band model
    approximation of the Exponential Wide Band model. The fixed-bands implementation
    pre-computed band limits for a set of temperature, pressure, species densities and
    path length representative of fire scenarios. The resulting bands are recorded in
    dictionary boxModelProperties, to be reused for predicting gas radiation properties
    in fires.

    Ref:
          "Radiative Heat Transfer", 2013, by Michael F. Modest
      Section 11.10, table 11.3

    box Model properties dictionary (must be in file: constant/multimediaRadProperties):

    allBands dictionary outlines all bands to be solved in fvDOM;
    all of these bands will be solved and properties will be required for each participating medium.
    When solving gas only, all but one of the transparent bands may be commented in allBands; the first
    band must be transparent (band limits here are not used; calculation lumps all transparent bands into
    a single composite, and energy fraction = 1 - Sum(energyFractions_of_radiatingBands)

    \verbatim
    allBands
    {
     bandLimits
     (
         (0.0 1.4e-6)
         (1.4e-6 2.49e-6)
         (2.49e-6 2.85e-6)
         (2.85e-6 4.15e-6)
         (4.15e-6 4.59e-6)
         (4.59e-6 5.49e-6)
         (5.49e-6 7.26e-6)
         (7.26e-6 20.0e-6)
     );
    }

    boxRadiatingSpecies
    {
        H2O
    {
        m          3;
        g_k        (1 1 1);
        eta_k      (3652 1595 3756);
        }
    CO2
    {
        m          3;
        g_k        (1 2 1);
        eta_k      (1351 666 2396);
        }
    }

    boxRadiatingBands
    {
        band1
    {
        bandLimits (1.3695e-6  1.3893e-6);
        species    H2O;
        eta_c      7250; // 1.38microns
        delta_k    (1 0 1);
        alpha_0    2.5;
        lambda     1.38;
        }
    band2
    {
            bandLimits (1.8481e-6  1.8907e-6);
        // overlap of one CO2 band at 5200 cm-1 and one H2O band at 5350 cm-1
        species    H2O;
        eta_c      5350; // 1.87microns
        delta_k    (0 1 1);
        alpha_0    3;
        lambda     1.87;
        overlap
        {
            species  CO2;
        alpha_0  0.06;
        eta_c    5200;
        lambda   2.0;
        }
        }
    band3
    {
        bandLimits (2.4907e-6  2.8531e-6);
        // overlap of one CO2 band at 3660 cm-1 and three H2O bands at 3760 cm-1
        species    H2O;
        eta_c      3760; // 2.7microns
        delta_k    (0 0 1);
        alpha_0    (23.4 2.3 0.2);
        lambda     2.7;
        overlap
        {
            species  CO2;
        alpha_0  4.0;
        eta_c    3660;
        lambda   2.7;
        }
        }
    band4
    {
        bandLimits (4.1494e-6  4.5851e-6);
        species    CO2;
        eta_c      2410; // 4.3microns
        delta_k    (0 0 1);
        alpha_0    110;
        lambda     4.3;
        }
    band5
        {
        bandLimits (5.4885e-6  7.2569e-6);
        species    H2O;
        eta_c      1600; // 6.3microns
        delta_k    (0 1 0);
        alpha_0    41.2;
        lambda     6.3;
        }
    band6
    {
        bandLimits (13.3156e-6  17.1821e-6);
        species    CO2;
        eta_c      667; // 15.0microns
        delta_k    (0 1 0);
        alpha_0    19;
        lambda     15;
        }

    }

    \endverbatim


SourceFiles
    boxFixedBandsAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef boxFixedBandsAbsorptionEmission_H
#define boxFixedBandsAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "absorptionCoeffs.H"
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class boxFixedBandsAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class boxFixedBandsAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //-  Modest's EWBM coefficients
        List<absorptionCoeffs> coeffs_;

        //- Bandwidths in units of wavenumbers (cm^-1)
        scalarList deltaEta_;

        //- Enum type for local or domain-averaged kappa, or domain-averaged density abs coeff, in box model
        enum absCoeffType
        {
            absCoeffAverage,
            absCoeffLocal,
            absCoeffAverageAlpha
        };

private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Name of file containing the box model dictionary
        const word boxFileName_;

        //- box model species dictionary
        const dictionary boxSpeciesDict_;

        //- box model properties dictionary
        const dictionary boxDict_;

        //- Thermo package
        const reactingMixture<gasHThermoPhysics>& gasThermo_;

        //- pressure
        const volScalarField& p_;

        //- temperature
        const volScalarField& T_;

        //- Species mass fractions
        const PtrList<volScalarField>& Y_;

        //- Total number of species in the mixture
        const label mixNspecies_;

        //- Species index of CO2 in mixture
        const label indexCO2_;

        //- Species index of H2O in mixture
        const label indexH2O_;

        //- Molecular weight of CO2
        const scalar molWtCO2_;

        //- Molecular weight of H2O
        const scalar molWtH2O_;

        //- Total number of radiating species
        label radNspecies_;

        //- eta_k properties of radiating species H2O and CO2
        scalarList etakH2O_;
        scalarList etakCO2_;

        //- Total number of bands
        label nBands_;

        //- Band limits
        List<Vector2D<scalar>> iBands_;

        //- Transparent/non-radiating bands
        labelList transparentBand_;

        //- Access to soot coefficient from case file
        const scalar Csoot_;

        //- Flag to indicate whether soot model is mixtureFractionSoot
        const Switch mixtureFractionSoot_;

        //- Get thermo refs, indexes and values
        const psiReactionThermo& thermoRef(const fvMesh& mesh) const;

        //- Get box model properties dictionary
        const dictionary boxDict(const fvMesh& mesh, const word key);

        //- Min. Temperature to include in domain averaging
        const scalar T_threshold_;

        //- Universal gas constant [J/(mol.K)]
        const scalar gasConstant_ = 8.3144598;

        //- Reference temperature in the Exponential Wide Band model
        const scalar T0_ = 100.;

        //- Planck second radiation constant ( hc/k_b [=] cm.K )
        const scalar C2_ = 1.438786;


protected:
    // Protected data

        //- Enumerator that identifies box model method as local or domain-average
        absCoeffType method_;


public:

        //- Runtime type information
        TypeName("boxFixedBandsAbsorptionEmission");

        // Constructors

            //- Construct from components
            boxFixedBandsAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


        //- Destructor
        virtual ~boxFixedBandsAbsorptionEmission(); // destructor

       // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> aCont(const label bandI = 0) const;

                    //- Auxiliary functions for computing aCont
                    inline scalar alphaFactor(const label eta_c, const word speciesName, const scalar T_avg) const;
                    inline scalar alpha(const scalar Temp, const label bandI = 0) const;
                    inline scalar rhoAlpha(const scalar rhoH2O, const scalar rhoCO2, const scalar alpha,
                    const scalar overlapAlpha, const label bandI = 0) const;

            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;

            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;

            // Gray model indicator
            inline bool isGrey() const
            {
                return false;
            }

            //- Number of bands
            inline label nBands() const
            {
                return nBands_;
            }

            //- Lower and upper limits of band i
            inline const Vector2D<scalar>& bands(const label i) const
            {
                return iBands_[i];
	    }

            //- Correct rays
            virtual void correct
	    (
		volScalarField& a_,
		PtrList<volScalarField>& aLambda
	     ) const;

            //- Correct energy fractions
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
