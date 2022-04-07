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
    Foam::thermoSingleLayerRP

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as Cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15 K

SourceFiles
    thermoSingleLayerRP.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayerRP_H
#define thermoSingleLayerRP_H

#include "thermoSingleLayer.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class heatTransferModel;
class phaseChangeModel;
class massAbsorptionModel;
class filmRadiationModel;

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayerRP Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayerRP
:
    public thermoSingleLayer
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayerRP(const thermoSingleLayerRP&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayerRP&);

        ////- Return boundary types for sensible enthalpy field
        //wordList hsBoundaryTypes();

protected:

    // Protected data

	//- Using Roll-Paper Model
	bool rollPaperModel_;

        //- Burnt page number
	volScalarField paper_;

        //- Solid surface temperature
	volScalarField Tpaper_;

        //- Surface temperature
	volScalarField TpaperSurf_;

	//- Thermally Thin zone (1:thin, 0:thick, between 0 and 1: partially thin)
	volScalarField pthin_;

	//- Debug heat flux
	volScalarField pflux_;

	//- Debug blocking
	volScalarField pblock_;

	//- Find pyrolysis front height (using pthin_)
	volScalarField pthinHeight_;
		
	//- Directional weighting parameters
	scalar minPyHeight_;
	scalar delDistance_;
	scalar scaleUp_;
	scalar scaleDown_;

	//- New method to calculate minimum distance to pyrolysis front
	volScalarField pyState_;    // 0: no neighbour; 1: positive neighbour; -1: negtive neighbour
	volScalarField pyRpmin_;    // minimum distance to delamination front
	volVectorField pyCpmin_;    // closest cell positon from corrent cell
	volScalarField pyRpminScaled_;    // minimum distance to delamination front
	volScalarField pyPage_;    // page number of the closest py front

    // Protected member functions

        //- Correct the film surface and wall temperatures
        virtual void updateSurfaceTemperatures();

public:

    //- Runtime type information
    TypeName("thermoSingleLayerRP");


    // Constructors

        //- Construct from components
        thermoSingleLayerRP
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType,
            const bool readFields = true
        );

    //- Destructor
    virtual ~thermoSingleLayerRP();

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "thermoSingleLayerRPI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
