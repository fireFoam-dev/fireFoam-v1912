/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "filmDiversionForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmDiversionForce, 0);
addToRunTimeSelectionTable(force, filmDiversionForce, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmDiversionForce::filmDiversionForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    force(film),
    coeff_(readScalar(dict.subDict("filmDiversionCoeffs").lookup("coeff")))
{
    Info << "film diversion coefficient: " << coeff_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

filmDiversionForce::~filmDiversionForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> filmDiversionForce::correct(volVectorField& U)
{
    const volScalarField& charFrac = filmModel_.charFrac();
    const volScalarField& deltaf = filmModel_.delta(); // kvm
    const dimensionedScalar delta0("d0", dimLength, 1e-4); // kvm

    // tmp<volScalarField> // kvm
    //     tRatio(deltaf/delta0); // kvm
    // tRatio->min(1.0); // kvm

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce/dimArea*dimVolume));

    //scalar coeff_ = 1.0;
    dimensionedScalar coeff(dimensionedScalar("coeff",dimMass/dimTime/dimTime,-coeff_));

    // TODO: use estimated surface temperature here instead of mean film temperature
    // tfvm.ref() += coeff*fvc::grad(charFrac)*tRatio; // kvm
    tfvm.ref() += coeff*fvc::grad(charFrac); // kvm

    if(filmModel_.time().outputTime()){ // kvm
        tmp<volVectorField> tGradCharFrac(fvc::grad(charFrac));
        tGradCharFrac->write();
    } // kvm

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
