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
    corrugatedCardboardFreeburnFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles.

    If my temperature is T1, neighbour is T2:

    T1 > T2: my side becomes fixedValue T2 bc, other side becomes fixedGradient.


    Example usage:
        myInterfacePatchName
        {
            type compressible::corrugatedCardboardFreeburn;
            neighbourFieldName  T;
            K                   lookup;
            KName               K;

            neighbourFieldRadiativeName Qr; // or none. Use on solid side
            fieldRadiativeName  Qr;         // fluid. Use on fluid side

            value       uniform 300;
        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: K : heat conduction at patch. Gets supplied how to lookup/calculate K:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'basicThermo' : use basicThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use basicSolidThermo K()
    - 'directionalSolidThermo' directionalK()

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

    Note: lags interface data so both sides use same data.
    - problem: schedule to calculate average would interfere
    with standard processor swaps.
    - so: updateCoeffs sets both to same Twall. Only need to do
    this for last outer iteration but don't have access to this.

SourceFiles
    corrugatedCardboardFreeburnFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef corrugatedCardboardFreeburnFvPatchScalarField_H
#define corrugatedCardboardFreeburnFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
//#include "thermoSingleLayer.H"
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
   Class corrugatedCardboardFreeburnFvPatchScalarField
   Declaration
\*---------------------------------------------------------------------------*/

class corrugatedCardboardFreeburnFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

//    typedef Foam::regionModels::surfaceFilmModels::thermoSingleLayer
//        filmModelType;

    typedef Foam::regionModels::pyrolysisModels::pyrolysisModel
        pyrolysisModelType;

private:

    // Private data

        ////- Name of film region
        //const word filmRegionName_;

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

        //- Convective Scaling Factor (as determined by Prateep's tests)
        const scalar convectiveScaling_;
        const scalar convectiveCoefficient_;

        ////- Minimum delta film to be consired dry
        //const scalar filmDeltaDry_;

        ////- Maximum delta film to be consired wet
        //const scalar filmDeltaWet_;

        //- Emissivity
        scalarField emissivity_;

	//- Bound total heat flux
	const bool qBound_;

	//- Maximum heat flux peak
	const scalar qMax_;

    //- Extra heat flux;
    //dimensionedScalar qExtra_;
    //const scalar qExtra_;

        //- how to obtain radiative flux
        enum operationMode
        {
            radFromNeighbour,
            radFromMe,
            noRad,
            unknown
        };
        //static const NamedEnum<operationMode, 4> operationModeNames;

        //- Previous iteration mode
        operationMode oldMode_;

        ////- Retrieve film model from the database
        //const filmModelType& filmModel() const;

        //- Retrieve pyrolysis model from the database
        const pyrolysisModelType& pyrModel() const;


public:

    //- Runtime type information
    TypeName("compressible::corrugatedCardboardFreeburn");


    // Constructors

        //- Construct from patch and internal field
        corrugatedCardboardFreeburnFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        corrugatedCardboardFreeburnFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        corrugatedCardboardFreeburnFvPatchScalarField
        (
            const
            corrugatedCardboardFreeburnFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new corrugatedCardboardFreeburnFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        corrugatedCardboardFreeburnFvPatchScalarField
        (
            const corrugatedCardboardFreeburnFvPatchScalarField&,
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
                new corrugatedCardboardFreeburnFvPatchScalarField
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
