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
    radiationCoupledBaseFF

Description
    Common functions to emissivity. It gets supplied from lookup into a
    dictionary or calculated by the solidThermo:

    - 'lookup' : lookup volScalarField with name
    - 'solidThermo' : use solidThermo emissivity()
    - 'pyroModel' : use emissivity supplied by Pyrolysis model

SourceFiles
    radiationCoupledBaseFF.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCoupledBaseFF_H
#define radiationCoupledBaseFF_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"
#include "pyrolysisModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class radiationCoupledBaseFF Declaration
\*---------------------------------------------------------------------------*/

class radiationCoupledBaseFF
{
public:

    //- Type of supplied emissivity
    enum emissivityMethodType
    {
        SOLIDRADIATION,
        PYROLYSISMODELE,
        LOOKUPE
    };

    //- Type of supplied absorptivity
    enum absorptivityMethodType
    {
        EMISSIVITY,
        PYROLYSISMODELA,
        LOOKUPA,
        SOLIDRADIATIONA
    };

    typedef Foam::regionModels::pyrolysisModels::pyrolysisModel
        pyroModelType;


private:

    // Private data

        static const NamedEnum<emissivityMethodType, 3>
            emissivityMethodTypeNames_;

        static const NamedEnum<absorptivityMethodType, 4>
            absorptivityMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;


protected:

    // Protected data

        //- How to get emissivity
        const emissivityMethodType method_;

        //- How to get absorptivity
        const absorptivityMethodType absMethod_;

         //- Emissivity
         // Cached locally when is read from dictionary (lookup mode)
        scalarField emissivity_;

         //- Absorptivity
         // Cached locally when is read from dictionary (lookup mode)
        scalarField absorptivity_;

        //- Retrieve pyrolysis model from the database
        const pyroModelType& pyroModel(const fvMesh&,const word&) const;

public:

    // Constructors

        //- Construct from patch, emissivity mode and emissivity
        radiationCoupledBaseFF
        (
            const fvPatch& patch,
            const word& calculationMethodE,
            const word& calculationMethodA,
            const scalarField& emissivity,
            const scalarField& absorptivity
        );

        //- Construct from patch and dictionary
        radiationCoupledBaseFF
        (
            const fvPatch& patch,
            const dictionary& dict
        );


    // Member functions

        //- Method to obtain emissivity
        word emissivityMethod() const
        {
            return emissivityMethodTypeNames_[method_];
        }

        //- Method to obtain absorptivity
        word absorptivityMethod() const
        {
            return absorptivityMethodTypeNames_[absMethod_];
        }


        //- Calculate corresponding emissivity field
        scalarField emissivity();

        //- Calculate corresponding absorptivity field
        scalarField absorptivity();

        //- Calculate transmissivity field of film covering boundary
        scalarField filmTransmissivity();

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
