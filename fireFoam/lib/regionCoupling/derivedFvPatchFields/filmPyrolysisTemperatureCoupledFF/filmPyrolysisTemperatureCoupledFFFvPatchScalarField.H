/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::
    compressible::
    filmPyrolysisTemperatureCoupledFFFvPatchScalarField

Description

    Coupled temperature BC for solid-film-gas interface

    This same BC must be called in both the solid(pyrolysis) and gas(primary) regions
    on the coupled patch

    The film region mediates the BC value indirectly, but absorbing the heat flux
    terms from the gas region, where the film is present. 

    The solid temperature update is determined by the net heat flux calculated

    The gas temperature update is fixed value, based on a blend of the film 
    and pyrolysis region and the wetness fraction.


    Both the solid and gas phase BCs need dictionaries that look like the below
    example

    Example usage:
       myInterfacePatchName
       {
           type            compressible::filmPyrolysisTemperatureCoupledFF;
           regionType      fluid;       // fluid or solid
           filmRegion      filmRegion;  // The name of the region
           pyrolysisRegion fuelRegion;  // The name of the region
           Tnbr            T;           // Name of T field in nbr region
           K               none;        // Name of kappa in this region
           neighbourFieldRadiativeName none; 
           neighbourFieldConvectiveName none;
           fieldRadiativeName qin;
           fieldConvectiveName none;
           value           $internalField;
       }
        
    Note: K : heat conduction at patch. Gets supplied how to lookup/calculate K:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'basicThermo' : use basicThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use basicSolidThermo K()
    - 'directionalSolidThermo' directionalK()

SourceFiles
    filmPyrolysisTemperatureCoupledFFFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef filmPyrolysisTemperatureCoupledFFFvPatchScalarField_H
#define filmPyrolysisTemperatureCoupledFFFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "thermoSingleLayer.H"
#include "pyrolysisModel.H"

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
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class filmPyrolysisTemperatureCoupledFFFvPatchScalarField
   Declaration
\*---------------------------------------------------------------------------*/

class filmPyrolysisTemperatureCoupledFFFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase   // only used for evaluation of kappa
{
public:

    typedef Foam::regionModels::surfaceFilmModels::thermoSingleLayer
        filmModelType;

    typedef Foam::regionModels::pyrolysisModels::pyrolysisModel
        pyrolysisModelType;

private:

    // Private data

        //- Name of film region
        const word filmRegionName_;

        //- Name of pyrolysis region
        const word pyrolysisRegionName_;

        //- Name of field on the neighbour region
        const word TnbrName_;

         //- Name of the radiative heat flux on the neighbout region
        const word neighbourFieldRadiativeName_;
        const word neighbourFieldConvectiveName_;

        //- Name of the radiative heat flux in the my region
        const word fieldRadiativeName_;
        const word fieldConvectiveName_;

        //- Name of thermal conductivity field
        const word KName_;

        //- Emissivity
        scalarField emissivity_;

        //- Extra heat flux;
        //dimensionedScalar qExtra_;
        const scalar qExtra_;

        //- how to obtain radiative flux
        enum operationMode
        {
            radFromNeighbour,
            radFromMe,
            noRad,
            unknown
        };
        //static const NamedEnum<operationMode, 4> operationModeNames;

        //- Retrieve film model from the database
        const filmModelType& filmModel() const;

        //- Retrieve pyrolysis model from the database
        const pyrolysisModelType& pyrModel() const;


public:

    //- Runtime type information
    TypeName("compressible::filmPyrolysisTemperatureCoupledFF");


    // Constructors

        //- Construct from patch and internal field
        filmPyrolysisTemperatureCoupledFFFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmPyrolysisTemperatureCoupledFFFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        filmPyrolysisTemperatureCoupledFFFvPatchScalarField
        (
            const
            filmPyrolysisTemperatureCoupledFFFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new filmPyrolysisTemperatureCoupledFFFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        filmPyrolysisTemperatureCoupledFFFvPatchScalarField
        (
            const filmPyrolysisTemperatureCoupledFFFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new filmPyrolysisTemperatureCoupledFFFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding Qr field
        const scalarField& Qr() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
