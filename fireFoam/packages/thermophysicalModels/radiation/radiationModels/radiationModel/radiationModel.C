/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"
#include "fvmSup.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radiationModel, 0);
        defineRunTimeSelectionTable(radiationModel, T);
        defineRunTimeSelectionTable(radiationModel, dictionary);
    }
}

const Foam::word Foam::radiation::radiationModel::externalRadHeatFieldName_ =
    "qrExt";

const Foam::word Foam::radiation::radiationModel::primaryFluxName_ =
    "qprimaryRad";

 const Foam::word Foam::radiation::radiationModel::relfectedFluxName_ =
    "qreflective";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiation::radiationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::radiation::radiationModel::initialise()
{
    if (radiation_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

        if (this->found("absorptionEmissionModel"))
        {
            absorptionEmission_.reset
            (
                absorptionEmissionModel::New(*this, mesh_).ptr()
            );
        }

        if (this->found("scatterModel"))
        {
            scatter_.reset(scatterModel::New(*this, mesh_).ptr());
        }

        if (this->found("sootModel"))
        {
            soot_.reset(sootModel::New(*this, mesh_).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::radiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
//<fmglobal>
    T4fac_
    (
	IOobject
	(
	    "T4fac",
	    time_.timeName(),
	    mesh_,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh_,
	dimensionedScalar("one",dimless,1.0)
    ),
//</fmglobal>
    radiation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
//<fmglobal>
    T4fac_
    (
	IOobject
	(
	    "T4fac",
	    time_.timeName(),
	    mesh_,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh_,
	dimensionedScalar("one",dimless,1.0)
    ),
//</fmglobal>
    radiation_(lookupOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    if (readOpt() == IOobject::NO_READ)
    {
        radiation_ = false;
    }

    initialise();
}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
//<fmglobal>
    T4fac_
    (
	IOobject
	(
	    "T4fac",
	    time_.timeName(),
	    mesh_,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh_,
	dimensionedScalar("one",dimless,1.0)
    ),    
//</fmglobal>
    radiation_(lookupOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::radiationModel::read()
{
    if (regIOobject::read())
    {
        readEntry("radiation", radiation_);
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }

    return false;
}


void Foam::radiation::radiationModel::correct()
{
    if (!radiation_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }

    if (!soot_.empty())
    {
        soot_->correct();
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_*T4fac_ - 4.0*he/Cpv) // luwi
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    tmp<volScalarField> rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp.ref()
      - fvm::Sp(Rp()*pow3(T)/rhoCp.ref(), T)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    volScalarField& T
) const
{
    return
    (
        Ru()
      - fvm::Sp(Rp()*pow3(T), T)
    );
}


const Foam::radiation::absorptionEmissionModel&
Foam::radiation::radiationModel::absorptionEmission() const
{
    if (!absorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation absorptionEmission model, but model is "
            << "not activate" << abort(FatalError);
    }

    return *absorptionEmission_;
}


const Foam::radiation::sootModel&
Foam::radiation::radiationModel::soot() const
{
    if (!soot_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not activate" << abort(FatalError);
    }

    return *soot_;
}

/*
const Foam::radiation::transmissivityModel&
Foam::radiation::radiationModel::transmissivity() const
{
    if (!transmissivity_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not activate" << abort(FatalError);
    }

    return *transmissivity_;
}
*/
//<fmglobal: radiative energy source term (explicit) as volScalarField>
Foam::tmp<Foam::volScalarField> Foam::radiation::radiationModel::Sh() const
{
    const volScalarField T2(T_*T_);
    tmp<volScalarField> tQdot_rad
	(
	    new volScalarField
	    (
		IOobject
		(
		    "tQdot_rad",
		    mesh_.time().timeName(),
		    mesh_,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE,
		    false
		    ),
		-Rp()*T2*T2*T4fac_
	    )
	);

    tQdot_rad.ref().primitiveFieldRef() += Ru();

    return tQdot_rad;
}

Foam::tmp<Foam::volScalarField> Foam::radiation::radiationModel::updateT4fac()
{
    tmp<volScalarField> tT4fac
    (
	new volScalarField
	(
	    IOobject
	    (
		"tT4fac",
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE,
		false
	    ),
	    mesh_,
	    dimensionedScalar("one", dimless, 1.0)
	)
    );

    return tT4fac;
}
//</fmglobal>


// ************************************************************************* //
