/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "standardRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    standardRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardRadiation::standardRadiation
(
     surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, film, dict),
    qinPrimary_
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
        dimensionedScalar(dimMass/pow3(dimTime), Zero),
        film.mappedPushedFieldPatchTypes<scalar>()
    ),
    qemPrimary_
    (
        IOobject
        (
            "qem", // same name as qem on primary region to enable mapping
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero),
        film.mappedPushedFieldPatchTypes<scalar>()
    ),
    beta_(coeffDict_.lookupOrDefault<scalar>("beta",0.9)),
    kappaBar_(coeffDict_.lookupOrDefault<scalar>("kappaBar",1.e4))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardRadiation::correct()
{
    // Transfer qin from primary region
    qinPrimary_.correctBoundaryConditions();
    qemPrimary_.correctBoundaryConditions();
}


tmp<volScalarField> standardRadiation::Shs()
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
            dimensionedScalar(dimMass/pow3(dimTime), Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& Shs = tShs.ref();
    const scalarField& qinP = qinPrimary_; // Received at the coupled boundary from primary region
    const scalarField& qemP = qemPrimary_; // Emitted from the coupled boundary into primary region
                                           // , through the film (can also be absorped 

    const scalarField& delta = filmModel_.delta();

    // Fraction of radiation absorbed by film; 
    // based on ref: "Experimental study on radiation attenuation by a water film"
    //               D. Brissinger, G. Parent, P. Boulet, JQSRT 145 (2014) 160-168
    scalarField absoptivity = beta_*(1.0 - exp(-kappaBar_*delta));

    // The film itself is also an emitter
    const scalar sigma = constant::physicoChemical::sigma.value();
    const scalarField& Tem = qinPrimary_.db().lookupObject<volScalarField>("Tsf");
    scalarField qemF = beta_*sigma*pow4(Tem);

    Shs = (qinP - qemP)*absoptivity - qemF;

    return tShs;
}

tmp<volScalarField> standardRadiation::ShsConst() const
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
            dimensionedScalar(dimMass/pow3(dimTime), Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& Shs = tShs.ref();
    const scalarField& qinP = qinPrimary_; // Received at the coupled boundary from primary region
    const scalarField& qemP = qemPrimary_; // Emitted from the coupled boundary into primary region
                                           // , through the film (can also be absorped 

    const scalarField& delta = filmModel_.delta();

    // Fraction of radiation absorbed by film; 
    // based on ref: "Experimental study on radiation attenuation by a water film"
    //               D. Brissinger, G. Parent, P. Boulet, JQSRT 145 (2014) 160-168
    scalarField absoptivity = beta_*(1.0 - exp(-kappaBar_*delta));

    // The film itself is also an emitter
    const scalar sigma = constant::physicoChemical::sigma.value();
    const scalarField& Tem = qinPrimary_.db().lookupObject<volScalarField>("Tsf");
    scalarField qemF = beta_*sigma*pow4(Tem);

    Shs = (qinP - qemP)*absoptivity - qemF;

    return tShs;
}

const tmp<volScalarField> standardRadiation::transmissivity() const
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
