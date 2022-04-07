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

Class
    Foam::reactingOneDimFF

Description
    Reacting, 1-D pyrolysis model

SourceFiles
    reactingOneDimFF.C

\*---------------------------------------------------------------------------*/

#ifndef reactingOneDimFF_H
#define reactingOneDimFF_H

#include "pyrolysisModel.H"
#include "basicSolidChemistryModel.H"
#include "radiationModel.H"
#include "thermoSingleLayer.H"

#define DEBUG(x) {                                              \
                std::streamsize p = std::cout.precision();              \
                std::ios::fmtflags myFlags;                             \
                myFlags = cout.flags();                                 \
                std::cout.precision(10);                                \
                std::cout.setf(std::ios::fixed,std::ios::floatfield);   \
                std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
                std::cout << "p" << Pstream::myProcNo();                \
                std::cout << " " << #x " = " << x << std::endl;         \
                std::cout.precision(p);                                 \
                std::cout.flags(myFlags);                               \
            }
#define TRACE(s) {                                              \
                std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
                std::cout << "p" << Pstream::myProcNo();                \
                std::cout << " " << #s << std::endl;                    \
                s;                                                      \
            }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{


/*---------------------------------------------------------------------------*\
                      Class reactingOneDimFF Declaration
\*---------------------------------------------------------------------------*/

class reactingOneDimFF
:
    public pyrolysisModel
{

public:

    typedef Foam::regionModels::surfaceFilmModels::thermoSingleLayer
        filmModelType;

private:

    // Private member functions

        //- No copy construct
        reactingOneDimFF(const reactingOneDimFF&) = delete;

        //- No copy assignment
        void operator=(const reactingOneDimFF&) = delete;

        //- Read model controls
        void readReactingOneDimControls();

        //- Found film region?
        bool filmRegionFound();

        //- Retrieve film model from the database
        const filmModelType& filmModel() const;

protected:



    // Protected data

        //- Reference to solid thermo
        autoPtr<solidReactionThermo> solidThermo_;

        //- Reference to the solid chemistry model
        autoPtr<basicSolidChemistryModel> solidChemistry_;

        //- Pointer to radiation model
        autoPtr<radiation::radiationModel> radiation_;


        // Reference to solid thermo properties

            //- Density [kg/m3]
            volScalarField rho_;

            //- List of solid components
            PtrList<volScalarField>& Ys_;

            // Non-const access to enthalpy
            volScalarField& h_;


        // Solution parameters

            //- Number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Maximum diffusivity
            scalar maxDiff_;

            //- Minimum delta for combustion
            scalar minimumDelta_;

            //- Minimum solid temperature
            scalar Tmin_;

            //- Maximum solid temperature
            scalar Tmax_;


        // Fields

            //- Total gas mass flux to the primary region [kg/m2/s]
            surfaceScalarField phiGas_;

            //- Sensible enthalpy gas flux [J/m2/s]
            volScalarField phiHsGas_;

            //- Heat release rate [J/s/m3]
            volScalarField chemistryQdot_;

            //- Emissivity at the boundary
            volScalarField emmBnd_;        

            //- Absorptivity at the boundary
            volScalarField absBnd_;        

            //- Transmissivity of water film on the boundary
            volScalarField filmTransmissivity_;


        // Source term fields

            //- Coupled region radiative heat flux [W/m2]
            //  Requires user to input mapping info for coupled patches
            //volScalarField qrCoupled_;

            //- In depth radiative heat flux [W/m2]
            volScalarField qr_;

            //- char fraction (used in film model) [-]
            volScalarField charFrac_;


        // Checks

            //- Cumulative lost mass of the condensed phase [kg]
            dimensionedScalar lostSolidMass_;

            //- Cumulative mass generation of the gas phase [kg]
            dimensionedScalar addedGasMass_;

            //- Total mass gas flux at the pyrolysing walls [kg/s]
            scalar totalGasMassFlux_;

            //- Total heat release rate [J/s]
            dimensionedScalar totalHeatRR_;


        // Options

            //- Add gas enthalpy source term
            bool  gasHSource_;

            //- Add in depth radiation source term
            bool  qrHSource_;

            //- Use chemistry solvers (ode or sequential)
            bool useChemistrySolvers_;

            //- Name of the film region (optional)
            word filmRegionName_;

            //- Found film region
            bool filmRegionFound_;

        // CharOxi 

            //- Total char volume fraction [-], only computed at patch
            //volScalarField Xchar_;

            //- Total char volume fraction [kg], only computed at patch
            //volScalarField mChar_;

            //- Total char volume fraction [kg], only computed at patch
            volScalarField mCharBurnt_;

            //- Heat release due to char oxidation [J/s/m3]
            volScalarField charOxiQdot_;

            //- Total O2 mass flux from the primary region [kg/s]
            //volScalarField phiO2_;

            //- Total CO2 mass flux to the primary region  [kg/s]
            //volScalarField phiCO2_;

            //- Temperature when char start oxidation
            scalar TcharStart_;

            //- Temperature when char oxidation rate reach maximum
            scalar TcharOxiMax_;
 
            //- Index of char species
            label charIndex_;

            //- Does the species list include char?
            bool charring_;
            
            
            
     // Protected member functions

        //- Read control parameters from dictionary
        bool read();

        //- Read control parameters from dict
        bool read(const dictionary& dict);

        //- Update submodels
        virtual void updateFields();

        //- Update emissivity and absorptivity at boundaries
        virtual void updateBndEmmAbs();

        //- Update transmissivity of water film on boundaries
        virtual void updateFilmTransmissivity();

        //- Update/move mesh based on change in mass
        void updateMesh(const scalarField& mass0);

        //- Update radiative flux in pyrolysis region
        void updateqr();

        //- Update enthalpy flux for pyrolysis gases
        void updatePhiGas();

        //- Mass check
        void calculateMassTransfer();

        //- Update char oxidation rate
        void updateCharOxi();


        // Equations

            //- Solve continuity equation
            void solveContinuity();

            //- Solve energy
            virtual void solveEnergy();

            //- Solve solid species mass conservation
            virtual void solveSpeciesMass();


public:

    //- Runtime type information
    TypeName("reactingOneDimFF");


    // Constructors

        //- Construct from type name and mesh
        reactingOneDimFF
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name, mesh and dictionary
        reactingOneDimFF
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );




    //- Destructor
    virtual ~reactingOneDimFF();


    // Member Functions

        // Access

            //- Fields

                //- Return const density [Kg/m3]
                const volScalarField& rho() const;

                //- Return const temperature [K]
                virtual const volScalarField& T() const;

                //- Return specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const;

                //- Return the region absorptivity [1/m]
                virtual tmp<volScalarField> kappaRad() const;

                //- Return the region thermal conductivity [W/m/k]
                virtual tmp<volScalarField> kappa() const;

                //- Return the total gas mass flux to primary region [kg/m2/s]
                virtual const surfaceScalarField& phiGas() const;

        // Char index
            label charIndex() const;

        // Solution parameters

            //- Return the number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Return max diffusivity allowed in the solid
            virtual scalar maxDiff() const;



        // Helper functions

            //- External hook to add mass to the primary region
            virtual scalar addMassSources
            (
                const label patchi,            // patchi on primary region
                const label facei              // facei of patchi
            );

            //- Mean diffusion number of the solid region
            virtual scalar solidRegionDiffNo() const;


       // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the pyrolysis equations
            virtual void evolveRegion();

            //- Post-evolve region
            virtual void postEvolveRegion();


       // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactingOneDimFFI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
