/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "boundaryRadiationProperties.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(boundaryRadiationProperties, 0);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::boundaryRadiationProperties
(
    const fvMesh& mesh
)
:
    MeshObject
    <
        fvMesh,
        Foam::GeometricMeshObject,
        boundaryRadiationProperties
    >(mesh),
    radBoundaryPropertiesPtrList_(mesh.boundary().size())
{
    IOobject boundaryIO
    (
        boundaryRadiationProperties::typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (boundaryIO.typeHeaderOk<IOdictionary>(true))
    {
        const radiationModel& radiation =
            mesh.lookupObject<radiationModel>
            (
                "radiationProperties"
            );

        // Model number of bands
        label nBands = radiation.nBands();

        const IOdictionary radiationDict(boundaryIO);

        forAll(mesh.boundary(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (radiationDict.isDict(pp.name()))
            {
                const dictionary& dict = radiationDict.subDict(pp.name());

                radBoundaryPropertiesPtrList_[patchi].reset
                (
                    boundaryRadiationPropertiesPatch::New(dict, pp)
                );
#if 1//luwi--when we use 'lookup', banding is handled by model specified in IDefault
                if (radBoundaryPropertiesPtrList_[patchi]->nBands()<0)
                {
                    continue;
                }
#endif
                if (nBands != radBoundaryPropertiesPtrList_[patchi]->nBands())
                {
                    FatalErrorInFunction
                        << "Radiation bands : " <<  nBands << nl
                        << "Bands on patch : " << patchi << " is "
                        << radBoundaryPropertiesPtrList_[patchi]->nBands()
                        << abort(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * *  //

Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::emissivity
(
    const label patchi,
    const label bandi,
    vectorField* incomingDirection,
    scalarField* T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->e
        (
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::scalar Foam::radiation::boundaryRadiationProperties::faceEmissivity
(
    const label patchi,
    const label faceI,
    const label bandi,
    vector incomingDirection,
    scalar T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->e
        (
            faceI,
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return Zero;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::absorptivity
(
    const label patchi,
    const label bandi,
    vectorField* incomingDirection,
    scalarField* T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->a
        (
            bandi,
            incomingDirection,
            T
        );
    }

     FatalErrorInFunction
         << "Patch : " << mesh().boundaryMesh()[patchi].name()
         << " is not found in the boundaryRadiationProperties. "
         << "Please add it"
         << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::scalar Foam::radiation::boundaryRadiationProperties::faceAbsorptivity
(
    const label patchi,
    const label faceI,
    const label bandi,
    vector incomingDirection,
    scalar T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->a
        (
            faceI,
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return Zero;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::transmissivity
(
    const label patchi,
    const label bandi,
    vectorField* incomingDirection,
    scalarField* T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->t
        (
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::scalar Foam::radiation::boundaryRadiationProperties::faceTransmissivity
(
    const label patchi,
    const label faceI,
    const label bandi,
    vector incomingDirection,
    scalar T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->t
        (
            faceI,
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return Zero;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::diffReflectivity
(
    const label patchi,
    const label bandi,
    vectorField* incomingDirection,
    scalarField* T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->rDiff
        (
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::scalar Foam::radiation::boundaryRadiationProperties::faceDiffReflectivity
(
    const label patchi,
    const label faceI,
    const label bandi,
    vector incomingDirection,
    scalar T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->rDiff
        (
            faceI,
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return Zero;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::specReflectivity
(
    const label patchi,
    const label bandi,
    vectorField* incomingDirection,
    scalarField* T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->rSpec
        (
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::scalar Foam::radiation::boundaryRadiationProperties::faceSpecReflectivity
(
    const label patchi,
    const label faceI,
    const label bandi,
    vector incomingDirection,
    scalar T
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->rSpec
        (
            faceI,
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return Zero;
}


// ************************************************************************* //
