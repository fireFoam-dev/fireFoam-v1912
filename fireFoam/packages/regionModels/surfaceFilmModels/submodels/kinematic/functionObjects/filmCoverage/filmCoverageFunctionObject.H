/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::filmCoverageFunctionObject

Description
    Region model function object to output the film coverage on a per-patch
    basis

SourceFiles
    filmCoverageFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef filmCoverageFunctionObject_H
#define filmCoverageFunctionObject_H

#include "regionModelFunctionObject.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regionModels
{

namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                 Class filmCoverageFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class filmCoverageFunctionObject
:
    public regionModelFunctionObject
{
protected:

    // Protected data

        //- Patch names
        wordReList patchNames_;
        

public:

    //- Runtime type information
    TypeName("filmCoverage");

    // Constructors

        //- Construct from dictionary
        filmCoverageFunctionObject(const dictionary& dict, regionModel& film);

        //- Construct copy
        filmCoverageFunctionObject(const filmCoverageFunctionObject& ppm);

        //- Construct and return a clone
        virtual autoPtr<regionModelFunctionObject> clone() const
        {
            return autoPtr<regionModelFunctionObject>
            (
                new filmCoverageFunctionObject(*this)
            );
        }


    //- Destructor
    virtual ~filmCoverageFunctionObject();


    // I-O

        //- write
        virtual void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
