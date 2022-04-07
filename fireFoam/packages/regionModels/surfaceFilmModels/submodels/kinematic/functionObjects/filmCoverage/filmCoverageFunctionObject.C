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

\*---------------------------------------------------------------------------*/

#include "filmCoverageFunctionObject.H"
#include "wordReList.H"
#include "polyBoundaryMesh.H"
#include "surfaceFilmModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

namespace surfaceFilmModels
{


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmCoverageFunctionObject, 0);
addToRunTimeSelectionTable
(
    regionModelFunctionObject,
    filmCoverageFunctionObject,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmCoverageFunctionObject::filmCoverageFunctionObject
(
    const dictionary& dict,
    regionModel& film
)
:
    regionModelFunctionObject(dict, film, typeName),
    patchNames_(dict.lookup("patches"))
{}


filmCoverageFunctionObject::filmCoverageFunctionObject
(
    const filmCoverageFunctionObject& fcfo
)
:
    regionModelFunctionObject(fcfo),
    patchNames_(fcfo.patchNames_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

filmCoverageFunctionObject::~filmCoverageFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void filmCoverageFunctionObject::write() const
{
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel& film =
        dynamic_cast<const regionModels::surfaceFilmModels::surfaceFilmRegionModel&>
        (
            regionModel_
        );

    const polyBoundaryMesh& pbm = film.regionMesh().boundaryMesh();
    labelHashSet patchIDs(pbm.patchSet(patchNames_));

    const scalarField& magSf = film.magSf();
    const scalarField& alpha =  film.alpha();

    if (patchIDs.size())
    {
        Info<< nl << type() << " output:" << nl;
    }

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchI = iter.key();

        const polyPatch& pp = pbm[patchI];

        const labelUList& faceCells = pp.faceCells();

        scalar sumMagSf = 0.0;
        scalar sumAlphaMagSf = 0.0;
        forAll(faceCells, i)
        {
            const label cellI = faceCells[i];

            sumMagSf += magSf[cellI];
            sumAlphaMagSf += alpha[cellI]*magSf[cellI];
        }

        reduce(sumMagSf, sumOp<scalar>());
        reduce(sumAlphaMagSf, sumOp<scalar>());

        Info<< "    Patch: " << pp.name() << nl
            << "        area     : " << sumMagSf << nl
            << "        coverage : " << sumAlphaMagSf/sumMagSf << nl;
    }

    if (patchIDs.size())
    {
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
