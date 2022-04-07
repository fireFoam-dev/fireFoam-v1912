/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "unitConversion.H"
#include "fvPatchField.H"
#include "meshWavePatchDistMethod.H"
#include "kinematicSingleLayer.H" // kvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactAngleForce, 0);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void contactAngleForce::initialise()
{
    const wordRes zeroForcePatches
    (
        coeffDict_.lookupOrDefault<wordRes>("zeroForcePatches", wordRes())
    );

    if (zeroForcePatches.size())
    {
        const polyBoundaryMesh& pbm = filmModel_.regionMesh().boundaryMesh();
        const scalar dLim = coeffDict_.get<scalar>("zeroForceDistance");

        Info<< "        Assigning zero contact force within " << dLim
            << " of patches:" << endl;

        labelHashSet patchIDs = pbm.patchSet(zeroForcePatches);

        for (const label patchi : patchIDs)
        {
            Info<< "            " << pbm[patchi].name() << endl;
        }

        // Temporary implementation until run-time selection covers this case
        patchDistMethods::meshWave dist(filmModel_.regionMesh(), patchIDs);
        volScalarField y
        (
            IOobject
            (
                "y",
                filmModel_.regionMesh().time().timeName(),
                filmModel_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            filmModel_.regionMesh(),
            dimensionedScalar("y", dimLength, GREAT)
        );
        dist.correct(y);

        mask_ = pos0(y - dimensionedScalar("dLim", dimLength, dLim));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    const word& typeName,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    force(typeName, film, dict),
    Ccf_(coeffDict_.get<scalar>("Ccf")),
    rndGen_(),
    mask_
    (
        IOobject
        (
            typeName + ":contactForceMask",
            filmModel_.time().timeName(),
            filmModel_.regionMesh(),
            // IOobject::MUST_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        filmModel_.regionMesh(),
        dimensionedScalar("mask", dimless, 1.0)
        // zeroGradientFvPatchScalarField::typeName
     ),
    nHits_ // kvm
    (
        IOobject
        (
            typeName + ":nHits",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

contactAngleForce::~contactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> contactAngleForce::correct(volVectorField& U)
{
    const kinematicSingleLayer& film =
        dynamic_cast<const kinematicSingleLayer&>(filmModel_); // kvm

    tmp<volVectorField> tForce
    (
        new volVectorField
        (
            IOobject
            (
                typeName + ":contactForce",
                filmModel_.time().timeName(),
                filmModel_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            filmModel_.regionMesh(),
            dimensionedVector(dimForce/dimArea, Zero)
        )
    );

    vectorField& force = tForce.ref().primitiveFieldRef();

    const labelUList& own = filmModel_.regionMesh().owner();
    const labelUList& nbr = filmModel_.regionMesh().neighbour();

    const scalarField& magSf = filmModel_.magSf();

    const volScalarField& alpha = filmModel_.alpha();
    const volScalarField& sigma = filmModel_.sigma();
    const volScalarField& deltaf = film.delta();

    volVectorField gradAlpha(fvc::grad(alpha));

    const tmp<volScalarField> ttheta = theta();
    const volScalarField& theta = ttheta();

    scalarField nHits(film.regionMesh().nCells(), 0.0);

    // minimum film thickness used for scaling contact angle force
    const scalar deltaf0 = 2e-4;

    // list of all faces, with own and nbr being the cell pair
    forAll(nbr, facei)
    {
        const label cellO = own[facei];
        const label cellN = nbr[facei];

        label celli = -1;
        label cellOther = -1;
        if ((alpha[cellO] > 0.5) && (alpha[cellN] < 0.5))
        {
            celli = cellO;
            cellOther = cellN;
        }
        else if ((alpha[cellO] < 0.5) && (alpha[cellN] > 0.5))
        {
            celli = cellN;
            cellOther = cellO;
        }

        // const scalar dxInverse = film.regionMesh().deltaCoeffs()[facei];
        // const scalar area = film.regionMesh().magSf()[facei];
        if (celli != -1)
        {
            // is this the right length to be using?
            const scalar dxInverse = film.regionMesh().deltaCoeffs()[facei]; // really 1/dxInverse
            // const scalar area = film.regionMesh().magSf()[facei]; // really 1/dxInverse
            const vector n =
                gradAlpha[celli]/(mag(gradAlpha[celli]) + ROOTVSMALL);
            // won't the contact angle be changing each timestep?
            // scalar theta = cos(degToRad(distribution_->sample()));
            scalar ratio = min(deltaf[celli]/deltaf0,1.0);
            scalar cosTheta = cos(theta[celli]);
            // force[celli] += mask_[celli]*Ccf_*n*sigma[celli]*(1.0 - cosTheta)/dxInverse;
            force[celli] += mask_[celli]*n*sigma[celli]*(1.0 - cosTheta)/Ccf_*ratio; // using Ccf_ as characteristic length
            // TODO: seems like we want to have the contact angel force on both sides of
            // the contact line interface.  This way, the residual liquid will be forced out.
            // (not sure about this one)
            // force[cellOther] += mask_[celli]*Ccf_*n*sigma[celli]*(1.0 - cosTheta)/dxInverse;
            // force[celli] += mask_[celli]*Ccf_*n*sigma[celli]*(1.0 - cosTheta)*dxInverse; //bug fix
            // force[celli] += mask_[celli]*Ccf_*n*sigma[celli]*(1.0 - cosTheta)*area; //bug fix
            nHits[celli]++;
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchField<scalar>& alphaf = alpha.boundaryField()[patchi];
        const scalarField& dxInverse = alphaf.patch().deltaCoeffs();
        // const scalarField& area = alphaf.patch().magSf();
        // why do we loop over the top and bottom patches?
        const labelUList& faceCells = alphaf.patch().faceCells();

        forAll(alphaf, facei)
        {
            label cellO = faceCells[facei];

            if ((alpha[cellO] > 0.5) && (alphaf[facei] < 0.5))
            {
                // won't the contact angle be changing each timestep?
                const vector n =
                    gradAlpha[cellO]/(mag(gradAlpha[cellO]) + ROOTVSMALL);
                // scalar cosTheta = cos(degToRad(distribution_->sample()));
                // ratio is important for code stability
                // when alpha=1 and deltaf is very small, the contact angle force
                // can be so large as to cause instabilities in the velocity.
                scalar ratio = min(deltaf[cellO]/deltaf0,1.0);
                scalar cosTheta = cos(theta[cellO]);
                // is this the right dxInverse to be using?
                // force[cellO] += mask_[cellO]*Ccf_*n*sigma[cellO]*(1.0 - cosTheta)/dxInverse[facei];
                force[cellO] += mask_[cellO]*n*sigma[cellO]*(1.0 - cosTheta)/Ccf_*ratio;
                // force[cellO] += mask_[cellO]*Ccf_*n*sigma[cellO]*(1.0 - cosTheta)*area[facei];
                nHits[cellO]++;
            }
        }
    }

    // why do we divide by nHits?  doesn't make sense
    // force /= (max(nHits, scalar(1.0))*magSf);
    // force /= magSf;
    tForce.ref().correctBoundaryConditions();

    if (filmModel_.regionMesh().time().writeTime())
    {
        tForce().write();
        nHits_.primitiveFieldRef() = nHits;
    }

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce/dimArea*dimVolume)); // why is this *dimVolume?

    tfvm.ref() += tForce; // doesn't this need to be *volume?

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
