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
    Foam::massAbsorptionModel

Description
    Base class for surface film phase change models

SourceFiles
    massAbsorptionModel.C
    massAbsorptionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef massAbsorptionModel_H
#define massAbsorptionModel_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class massAbsorptionModel Declaration
\*---------------------------------------------------------------------------*/

class massAbsorptionModel
:
    public filmSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        massAbsorptionModel(const massAbsorptionModel&);

        //- Disallow default bitwise assignment
        void operator=(const massAbsorptionModel&);


protected:

    // Protected Member Functions

        //- Latest mass transfer due to phase change
        scalar latestMassAbs_;

        //- Total mass transfer due to phase change
        scalar totalMassAbs_;


public:

    //- Runtime type information
    TypeName("massAbsorptionModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             massAbsorptionModel,
             dictionary,
             (
                surfaceFilmRegionModel& owner,
                const dictionary& dict
             ),
             (owner, dict)
         );

    // Constructors

        //- Construct null
        massAbsorptionModel(surfaceFilmRegionModel& owner);

        //- Construct from type name, dictionary and surface film model
        massAbsorptionModel
        (
            const word& modelType,
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<massAbsorptionModel> New
        (
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~massAbsorptionModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const scalar dt,
                scalarField& availableMass,
                volScalarField& dMass,
                volScalarField& dEnergy
            );

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            ) = 0;


        // I-O

            //- Provide some feedback
            virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
