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
    Foam::thermoSingleLayerRollPaper

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as Cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15 K

SourceFiles
    thermoSingleLayerRollPaper.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayerRollPaper_H
#define thermoSingleLayerRollPaper_H

#include "thermoSingleLayer.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayerRollPaper Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayerRollPaper
:
    public thermoSingleLayer
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayerRollPaper(const thermoSingleLayerRollPaper&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayerRollPaper&);

        ////- Return boundary types for sensible enthalpy field
        //wordList hsBoundaryTypes();

protected:

    // Protected data

        //- Burnt page number
	//volScalarField paper_;
	volScalarField filmPaperID_;

        //- Solid surface temperature
	//volScalarField Tpaper_;
	volScalarField filmTpaper_;

        //- Surface temperature
	//volScalarField TpaperSurf_;
	volScalarField filmTsurface_;

	//- Thermally Thin zone (1:thin, 0:thick, between 0 and 1: partially thin)
	//volScalarField pthin_;
	volScalarField filmPaperPeelingZone_;

	//- Extra heat flux on first paper at pyrolysis front
	//volScalarField pHeatFluxExtra_;

	//- Heat flux ratio (corrected heat flux / original heat flux)
	//volScalarField pHeatFluxRatio_;

	//- Debug heat flux
	//volScalarField pflux_;
	volScalarField filmQnet_;

	//- Debug blocking
	//volScalarField pblock_;
	volScalarField filmBlockFactor_;

	//- Debug mass loss rate
	//volScalarField pMLR_;
	volScalarField filmMassReleaseRate_;

	//- Find pyrolysis front height (using pthin_)
	//volScalarField pthinHeight_;
	volScalarField filmPaperPyrolysisFront_;

	//- Number of neighbor cells with positive pageID difference
	//volScalarField pDelPositive_;
	volScalarField filmPositiveDelNeighbor_;
		
	//- Number of neighbor cells with negtive pageID difference
	volScalarField filmNegtiveDelNeighbor_;
		
	//- Mimimum pyrolysis front height to activate peeling model
	scalar filmMinPyrolysisHeight_;

	//- Peeling distance at pyrolysis front
	scalar filmPeelingDistance_;

	//- Distance to closest pyrolysis front
	volScalarField filmMinPyrolysisDistance_;

	//- Cell position on pyrolysis front with the minimum distance
	volVectorField filmPyrolysisFrontCellPosition_;

	//- Page ID of the corresponding pyrolysis front
	volScalarField filmPyrolysisFrontPaperID_;

	//- Turn edge effect when reach a certain number of burnt paper
	scalar paperStackEffectID_;

    // Protected member functions

        //- Correct the film surface and wall temperatures
        virtual void updateSurfaceTemperatures();

	//- Correct roll paper thermally thin zone
	void updateRollPaperThermallyThinZone();

public:

    //- Runtime type information
    TypeName("thermoSingleLayerRollPaper");


    // Constructors

        //- Construct from components
        thermoSingleLayerRollPaper
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );

    //- Destructor
    virtual ~thermoSingleLayerRollPaper();

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "thermoSingleLayerRollPaperI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
