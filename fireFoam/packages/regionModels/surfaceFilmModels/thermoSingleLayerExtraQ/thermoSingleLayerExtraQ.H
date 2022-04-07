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
    Foam::thermoSingleLayerExtraQ

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as Cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15 K

SourceFiles
    thermoSingleLayerExtraQ.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayerExtraQ_H
#define thermoSingleLayerExtraQ_H

#include "thermoSingleLayerFmBox.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayerExtraQ Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayerExtraQ
:
    public thermoSingleLayerFmBox
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayerExtraQ(const thermoSingleLayerExtraQ&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayerExtraQ&);

        ////- Return boundary types for sensible enthalpy field
        //wordList hsBoundaryTypes();

protected:

    // Protected data

        //- Mass loss rate in pyrolysis zone (kg/m2/s)
	    volScalarField massLossRate_;

        //- Heat flux from gas-phase (qr + convectiveHeatFlux_T)
	    volScalarField totalHeatFlux_;

        //- Solid surface temperature
	    //volScalarField TsolidSurface_;

        //- Heat flux correction zone
	    volScalarField correctionZone_;
	    volScalarField burntZone_;

        //- Value from gas phase
        //surfaceScalarField phiSolid_;
        //surfaceScalarField qConvect_;
        volScalarField phiSolid_;
        volScalarField qConvect_;
        volScalarField qRadR_;

        //- Water absorption
        volScalarField waterAbsorbed_; // mm
        scalar absorptionRate_;        // mm/s
        scalar saturationThickness_;   // mm

        //scalar preHeatingDistance_;
        scalar preHeatingCells_;
        scalar preHeatingZoneMaxHeight_;

    // Protected member functions

        //- Correct the film surface and wall temperatures
        virtual void updateSubmodels();

	    //- Heat flux correction zone (pyrolysis front)
	    void getHeatFluxCorrectionZone();

public:

    //- Runtime type information
    TypeName("thermoSingleLayerExtraQ");


    // Constructors

        //- Construct from components
        thermoSingleLayerExtraQ
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );

    //- Destructor
    virtual ~thermoSingleLayerExtraQ();

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "thermoSingleLayerExtraQI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
