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
    Foam::thermoCloud

Description
    Virtual abstract base class for templated ThermoCloud

SourceFiles
    thermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef thermoCloud_H
#define thermoCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class thermoCloud Declaration
\*---------------------------------------------------------------------------*/

class thermoCloud
{
public:

    //- Runtime type information
    TypeName("thermoCloud");


    // Constructors

        //- Null constructor
        thermoCloud() = default;


    //- Destructor
    virtual ~thermoCloud() = default;


    // Member Functions

        // Access

            // Radiation

                //- Return tmp equivalent particulate emission
                virtual tmp<volScalarField> Ep(const label bandI) const = 0; // ankur

                //- Return tmp equivalent particulate absorption
                virtual tmp<volScalarField> ap(const label bandI) const = 0; // ankur

                //- Return tmp equivalent particulate scattering factor
                virtual tmp<volScalarField> sigmap(const label bandI) const = 0; // ankur

                //- Return tmp phase-function 
                virtual tmp<volScalarField> phaseFunc(const label bandI, const label sour, const label dest) const = 0; // ankur

                //- Return tmp asymmetry factor
                virtual tmp<volScalarField> asyFac(const label bandI) const = 0; // luwi

                //- Return tmp phase-function: Henyey-Greenstein computed directly each time 
                virtual tmp<volScalarField> phaseFuncHGDirect(const label bandI, const label sour, const label dest) const = 0; // luwi

                //- Return flag to indicate whether radiation is "on" for cloud
                virtual bool radiation() const = 0; // luwi

            // Sprinkler skipping
            
                //- Return the particle temperature scalar field
                //  Note: for particles belonging to this cloud only
                virtual const tmp<volScalarField> sprayTemperature() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
