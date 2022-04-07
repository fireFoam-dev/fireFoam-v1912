/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

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

#include "constantRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    constantRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantRadiation::constantRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, film, dict),
    qrConst_
    (
        IOobject
        (
            typeName + ":qrConst",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    ),
    mask_
    (
        IOobject
        (
            typeName + ":mask",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    qin_ // kvm
    (
        IOobject
        (
            "qin", // same name as qin on primary region to enable mapping
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
        zeroGradientFvPatchScalarField::typeName
        ),
    timeStart_(readScalar(coeffDict_.lookup("timeStart"))),
    duration_(readScalar(coeffDict_.lookup("duration"))),
    beta_(coeffDict_.lookupOrDefault<scalar>("beta",0.9)),
    kappaBar_(coeffDict_.lookupOrDefault<scalar>("kappaBar",1.e4))
{
    mask_ = pos0(mask_ - 0.5); // kvm
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantRadiation::~constantRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void constantRadiation::correct()
{
}


tmp<volScalarField> constantRadiation::Shs()
{
    tmp<volScalarField> tShs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Shs",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    const scalar time = film().time().value();

    const scalarField& Tem = qin_.db().lookupObject<volScalarField>("Tsf"); 

    const scalarField& delta = filmModel_.delta();

    if ((time >= timeStart_) && (time <= timeStart_ + duration_))
    {
        scalarField& Shs = tShs.ref();
        const scalarField& qr = qrConst_;
        // const scalarField& alpha = film.alpha();

        const scalar sigma = constant::physicoChemical::sigma.value();

        scalarField absorptivity = beta_*(1.0 - exp(-kappaBar_*delta));

        Info << "Radiation in: " << mask_[0]*qr[0]*absorptivity[0] << endl;
 
        Info << "Radiation out: " << sigma*absorptivity[0]*pow4(Tem[0]) << endl;

        // Shs = mask_*qr*alpha*absorptivity_;
        Shs = qr*mask_*absorptivity - beta_*sigma*pow4(Tem);
        //Shs = mask_*qr*absorptivity_ - sigma*absorptivity_*pow4(Tem);

        Info << "Rad net: " << Shs[0] << endl;

        // qin_ is used in turbulentTemperatureRadiationCoupledMixedST 
        // boundary condition
        qin_.primitiveFieldRef() = mask_*qr;
        qin_.correctBoundaryConditions();
    }

    return tShs;
}

tmp<volScalarField> constantRadiation::ShsConst() const
{
    tmp<volScalarField> tShs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Shs",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    const scalar time = film().time().value();

    const scalarField& Tem = qin_.db().lookupObject<volScalarField>("Tsf"); 

    const scalarField& delta = filmModel_.delta();

    if ((time >= timeStart_) && (time <= timeStart_ + duration_))
    {
        scalarField& Shs = tShs.ref();
        const scalarField& qr = qrConst_;

        const scalar sigma = constant::physicoChemical::sigma.value();

        scalarField absorptivity = beta_*(1.0 - exp(-kappaBar_*delta));

        Shs = (qr*mask_ - sigma*pow4(Tem))*absorptivity;
    }

    return tShs;
}


const tmp<volScalarField> constantRadiation::transmissivity() const
{
    tmp<volScalarField> ttau
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "transmissivity",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

	const scalarField& delta = filmModel_.delta();
    scalarField absorptivity = beta_*(1.0 - exp(-kappaBar_*delta));
	scalarField& tau = ttau.ref();

	tau = max(1.0 - absorptivity, SMALL); 

    return ttau;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
