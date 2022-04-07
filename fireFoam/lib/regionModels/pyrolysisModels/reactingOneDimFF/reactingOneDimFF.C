/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "reactingOneDimFF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvcLaplacian.H"
#include "absorptionEmissionModel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactingOneDimFF, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDimFF, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDimFF, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void reactingOneDimFF::updateCharOxi()
{


    scalar rhoChar = solidThermo_->composition().rho(charIndex(), 10000, 298);

    //Info << "mw of char is " <<  solidThermo_->composition().W(charIndex);

    // hardcode molecular weight for C and O2
    scalar mWO2 = 32;
    scalar mWChar = 12;
    scalar mWCO2 = 44;

    // hardcode HoC for char [J/kg]
    dimensionedScalar HocChar("HocChar",dimEnergy/dimMass,32.8e6);

    label localPyrolysisFaceI = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const fvPatch& patch = regionMesh().boundary()[patchi];

        scalarField Xcharp(patch.size(),Zero); 
        scalarField mCharp(patch.size(),Zero);
        scalarField phiO2p(patch.size(),Zero); 
        scalarField phiCO2p(patch.size(),Zero);

        scalarField& mCharBurntp = mCharBurnt_.boundaryFieldRef()[patchi];
        scalarField& charOxiQdotp = charOxiQdot_.boundaryFieldRef()[patchi];

        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
            patch.patch()
        );
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[mpp.samplePolyPatch().index()];

        fvPatchScalarField O2 = nbrPatch.lookupPatchField<volScalarField, scalar>("O2");
        scalarField alpha = nbrPatch.lookupPatchField<volScalarField, scalar>("thermo:alpha");    //for now use molecular alpha
        scalarField O2Int(O2.patchInternalField());
        scalarField alphaDelta(alpha * nbrPatch.deltaCoeffs());

        mpp.distribute(O2);
        mpp.distribute(O2Int);
        mpp.distribute(alphaDelta);

        phiO2p = - alphaDelta * (O2Int - 0.0) * patch.magSf(); //[kg/s] negative
        phiCO2p = - phiO2p * mWCO2 / mWO2;  // positive, same as phiGas

        scalarField deltaMO2(- phiO2p * time_.deltaT().value());   //[kg]

        const scalarField& cellV = regionMesh().V();

        forAll(Xcharp, faceI)
        {
            const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];
            scalar charVolume = 0.0;
            scalar hotCharVolume = 0.0;
            scalar totalVolume = 0.0;
            forAll(cells, k)
            {
                const label cellI = cells[k];
                scalar X = Ys_[charIndex()][cellI] * rho_[cellI] / rhoChar;

                charVolume += X * cellV[cellI];

                scalar hotCharFraction(max(0,min(1,(T()[cellI]-TcharStart_)/(TcharOxiMax_-TcharStart_))));

                //if(T()[cellI] > TcharOxidation_)    //- Turn off char oxidation if temperature too low
                //{
                //    hotCharVolume += X * cellV[cellI];
                //}

                hotCharVolume += hotCharFraction * X * cellV[cellI];
                totalVolume += cellV[cellI];
            }
            Xcharp[faceI] = charVolume / totalVolume;
            //mCharp[faceI] = rhoChar * charVolume;
            mCharp[faceI] = rhoChar * hotCharVolume;

            scalar charAvai = max(0, mCharp[faceI] - mCharBurntp[faceI]);
            scalar dmCharBurnt = 0.0;
            if (charAvai < deltaMO2[faceI] / mWO2 * mWChar)
            {
                dmCharBurnt = charAvai;
                phiO2p[faceI] = - dmCharBurnt / mWChar * mWO2;
                phiCO2p[faceI] = dmCharBurnt / mWChar * mWCO2;
            }
            else
            {
                dmCharBurnt = deltaMO2[faceI] / mWO2 * mWChar;
            }
            mCharBurntp[faceI] += dmCharBurnt;
            //Info<<"dbg-char: "<<time_.value()<<tab<<T()[cells[0]]<<tab<<totalVolume<<tab<<hotCharVolume<<tab
            //    <<Xcharp[faceI]<<tab<<mCharp[faceI]<<tab<<charAvai<<tab<<mCharBurntp[faceI]<<endl;

            //compute energy release rate in the first cell [kW/m3]
            charOxiQdot_[cells[0]] = HocChar.value() * dmCharBurnt
                                   / cellV[cells[0]]
                                   / time_.deltaT().value();


            //for HRR diagnostics (sum at the patch is the HRR) [kW]
            charOxiQdotp[faceI] = charOxiQdot_[cells[0]] * cellV[cells[0]];
            //Info << "charOxiQdot_" << charOxiQdot_[cells[0]] << endl;

            localPyrolysisFaceI++;
        }
    }
}

bool
reactingOneDimFF::filmRegionFound()
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return true;
        }
    }

    Info<< "NOTE: In reactingOneDimFF...\n"
        << "    No film region found for radiative absorption\n"
        << "    \"" << filmRegionName_
        << "\" will be treated as transparent"
        << nl;

    return false;
}

const reactingOneDimFF::filmModelType&
reactingOneDimFF::
filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }

    return **models.begin();
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void reactingOneDimFF::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.readEntry("nNonOrthCorr", nNonOrthCorr_);
    time().controlDict().readEntry("maxDi", maxDiff_);
    coeffs().readEntry("minimumDelta", minimumDelta_);
    Tmin_ = coeffs().getOrDefault("Tmin", 200.);
    Tmax_ = coeffs().getOrDefault("Tmax", 5000.);
    gasHSource_ = coeffs().getOrDefault("gasHSource", false);
    coeffs().readEntry("qrHSource", qrHSource_);
    useChemistrySolvers_ =
        coeffs().getOrDefault("useChemistrySolvers", true);
    filmRegionName_ =
        coeffs().getOrDefault<word>("filmRegion", "filmRegion");
    filmRegionFound_ = filmRegionFound();
}


bool reactingOneDimFF::read()
{
    if (pyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }

    return false;
}


bool reactingOneDimFF::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }

    return false;
}


void reactingOneDimFF::updateqr()
{
    // Update local qr from coupled qr field
    qr_ == dimensionedScalar(qr_.dimensions(), Zero);

    // Retrieve field from coupled region using mapped boundary conditions
    qr_.correctBoundaryConditions();

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        // qr is positive going in the solid
        // If the surface is emitting the radiative flux is set to zero
        qrBf[patchi] = max(qrBf[patchi], scalar(0));
    }

    const vectorField& cellC = regionMesh().cellCentres();

    tmp<volScalarField> kappa = kappaRad();

    // Propagate qr through 1-D regions
    label localPyrolysisFacei = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        const scalarField& qrp = qr_.boundaryField()[patchi];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchi];

        forAll(qrp, facei)
        {
            const scalar qr0 = qrp[facei];
            point Cf0 = Cf[facei];
            const labelList& cells = boundaryFaceCells_[localPyrolysisFacei++];
            scalar kappaInt = 0.0;
            forAll(cells, k)
            {
                const label celli = cells[k];
                const point& Cf1 = cellC[celli];
                const scalar delta = mag(Cf1 - Cf0);
                kappaInt += kappa()[celli]*delta;
                qr_[celli] = qr0*exp(-kappaInt);
                Cf0 = Cf1;
            }
        }
    }
}


void reactingOneDimFF::updatePhiGas()
{
    phiHsGas_ == dimensionedScalar(phiHsGas_.dimensions(), Zero);
    phiGas_ == dimensionedScalar(phiGas_.dimensions(), Zero);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas =
            solidChemistry_->gasHs(solidThermo_->p(), solidThermo_->T(), gasI);

        const volScalarField& HsiGas = tHsiGas();

        const volScalarField::Internal& RRiGas = solidChemistry_->RRg(gasI);

        surfaceScalarField::Boundary& phiGasBf = phiGas_.boundaryFieldRef();

        label totalFaceId = 0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchi = intCoupledPatchIDs_[i];

            scalarField& phiGasp = phiGasBf[patchi];
            const scalarField& cellVol = regionMesh().V();

            forAll(phiGasp, facei)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId++];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label celli = cells[k];
                    massInt += RRiGas[celli]*cellVol[celli];
                    phiHsGas_[celli] += massInt*HsiGas[celli];
                }

                phiGasp[facei] += massInt;

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchi
                        << " mass produced at face(local) : "
                        <<  facei
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
            }
        }
    }
}


void reactingOneDimFF::updateFields()
{
    if (qrHSource_)
    {
        updateqr();
    }

    //Note: Commented out (?) as the sensible gas energy is included in energy eq.
    updatePhiGas();

    if (charring_)
    {
        forAll(charFrac_,celli)
        {
            charFrac_[celli] = Ys_[charIndex()][celli];
        }
        charFrac_.correctBoundaryConditions();

        updateCharOxi();
    }
}


void reactingOneDimFF::updateBndEmmAbs()
{
    emmBnd_ == radiation_->absorptionEmission().e()();
    absBnd_ == radiation_->absorptionEmission().a()();
}

void reactingOneDimFF::updateFilmTransmissivity()
{
	const filmModelType& film = filmModel();
	forAll(intCoupledPatchIDs_, i)
	{
		const label patchI = intCoupledPatchIDs_[i];
        filmTransmissivity_.boundaryFieldRef()[patchI] =
			mapRegionPatchField<scalar>
            (
                film,
                "filmTransmissivity",
                patchI,
                true
            );
	}
}

void reactingOneDimFF::updateMesh(const scalarField& deltaV)
{
    if (!moveMesh_)
    {
        return;
    }

    Info<< "Initial/final volumes = " << gSum(regionMesh().V()) << ", "
        << gSum(regionMesh().V() + deltaV) << " [m3]" << endl;

    // Move the mesh
    const labelList moveMap = moveMesh(deltaV, minimumDelta_);

    // Flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 1)
        {
            solidChemistry_->setCellReacting(i, false);
        }
    }
}


void reactingOneDimFF::solveContinuity()
{
    DebugInFunction << endl;

    if (!moveMesh_)
    {
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho_) == -solidChemistry_->RRg()
        );

        rhoEqn.solve();
    }
    else
    {
        const scalarField deltaV
        (
            -solidChemistry_->RRg()*regionMesh().V()*time_.deltaT()/rho_
        );

        updateMesh(deltaV);
    }
}


void reactingOneDimFF::solveSpeciesMass()
{
    DebugInFunction << endl;

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi) == solidChemistry_->RRs(i)
        );

        if (regionMesh().moving())
        {
            surfaceScalarField phiYiRhoMesh
            (
                fvc::interpolate(Yi*rho_)*regionMesh().phi()
            );

            YiEqn -= fvc::div(phiYiRhoMesh);

        }

        YiEqn.solve(regionMesh().solver("Yi"));
        Yi.max(0.0);
        Yt += Yi;
    }

    Ys_[Ys_.size() - 1] = 1.0 - Yt;

}


void reactingOneDimFF::solveEnergy()
{
    DebugInFunction << endl;

    tmp<volScalarField> alpha(solidThermo_->alpha());

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_, h_)
      - fvm::laplacian(alpha, h_)
      + fvc::laplacian(alpha, h_)
      - fvc::laplacian(kappa(), T())
     ==
        chemistryQdot_
      + solidChemistry_->RRsHs()
      + charOxiQdot_
    );

    if (gasHSource_)
    {
        const surfaceScalarField phiGas(fvc::interpolate(phiHsGas_));
        hEqn += fvc::div(phiGas);
    }

    if (qrHSource_)
    {
        const surfaceScalarField phiqr(fvc::interpolate(qr_)*nMagSf());
        hEqn += fvc::div(phiqr);
    }

/*
    NOTE: The moving mesh option is only correct for reaction such as
    Solid -> Gas, thus the ddt term is compensated exactly by chemistrySh and
    the mesh flux is not necessary.
*/

    if (regionMesh().moving())
    {
        surfaceScalarField phihMesh
        (
            fvc::interpolate(rho_*h_)*regionMesh().phi()
        );

        hEqn -= fvc::div(phihMesh);
    }

    hEqn.relax();
    hEqn.solve();
}


void reactingOneDimFF::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchi]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistryQdot_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDimFF::reactingOneDimFF
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType),
    solidThermo_(solidReactionThermo::New(regionMesh())),
    solidChemistry_(basicSolidChemistryModel::New(solidThermo_())),
    radiation_(radiation::radiationModel::New(solidThermo_->T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_->rho()
    ),
    Ys_(solidThermo_->composition().Y()),
    h_(solidThermo_->he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),
    Tmin_(200),
    Tmax_(5000),
    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime, Zero)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime, Zero)
    ),

    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, Zero)
    ),
  
    emmBnd_
    (
        IOobject
        (
            "emmBnd",
            time().timeName(),
            regionMesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        radiation_->absorptionEmission().e()()
        //regionMesh(),
        //dimensionedScalar("zero", dimless/dimLength, 0.0)
        //dimensionedScalar("zero", dimless, 0.0)
    ),

    absBnd_
    (
        IOobject
        (
            "absBnd",
            time().timeName(),
            regionMesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        emmBnd_
        //regionMesh(),
        //dimensionedScalar("zero", dimless/dimLength, 0.0)
        //dimensionedScalar("zero", dimless, 0.0)
    ),

    filmTransmissivity_
    (
        IOobject
        (
            "filmTransmissivity",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchVectorField::typeName
    ),

    qr_
    (
        IOobject
        (
            "qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
        //dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        //zeroGradientFvPatchVectorField::typeName
    ),
    charFrac_
    (
        IOobject
        (
            "charFrac",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchVectorField::typeName
    ),
    lostSolidMass_(dimensionedScalar(dimMass, Zero)),
    addedGasMass_(dimensionedScalar(dimMass, Zero)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar(dimEnergy/dimTime, Zero)),
    gasHSource_(false),
    qrHSource_(false),
    useChemistrySolvers_(true),
    filmRegionName_("filmRegion"),
    filmRegionFound_(false),
    mCharBurnt_
    (
        IOobject
        (
            "mCharBurnt",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    charOxiQdot_
    (
        IOobject
        (
            "charOxiQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    TcharStart_(coeffs().lookupOrDefault<scalar>("TcharStart", 700)),
    TcharOxiMax_(coeffs().lookupOrDefault<scalar>("TcharOxiMax", 800)),
    charIndex_(-1),
    charring_(true)  // Will be set to false if 'char' not found in species list
{
    if (active_)
    {
        read();
    }

    charIndex_ = solidThermo_->composition().species()["char"];

    if (charIndex_ < 0)
    {
        charring_ = false;
    }

}


reactingOneDimFF::reactingOneDimFF
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType),
    solidThermo_(solidReactionThermo::New(regionMesh())),
    solidChemistry_(basicSolidChemistryModel::New(solidThermo_())),
    radiation_(radiation::radiationModel::New(solidThermo_->T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_->rho()
    ),
    Ys_(solidThermo_->composition().Y()),
    h_(solidThermo_->he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),
    Tmin_(200),
    Tmax_(5000),
    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime, Zero)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime, Zero)
    ),

    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, Zero)
    ),
    
    emmBnd_
    (
        IOobject
        (
            "emmBnd",
            time().timeName(),
            regionMesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        radiation_->absorptionEmission().e()()
        //regionMesh(),
        //dimensionedScalar("zero", dimless/dimLength, 0.0)
        //dimensionedScalar("zero", dimless, 0.0)
    ),

    absBnd_
    (
        IOobject
        (
            "absBnd",
            time().timeName(),
            regionMesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        emmBnd_
        //regionMesh(),
        //dimensionedScalar("zero", dimless/dimLength, 0.0)
        //dimensionedScalar("zero", dimless, 0.0)
    ),

    filmTransmissivity_
    (
        IOobject
        (
            "filmTransmissivity",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchVectorField::typeName
    ),

    qr_
    (
        IOobject
        (
            "qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    charFrac_
    (
        IOobject
        (
            "charFrac",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchVectorField::typeName
    ),
    lostSolidMass_(dimensionedScalar(dimMass, Zero)),
    addedGasMass_(dimensionedScalar(dimMass, Zero)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar(dimEnergy/dimTime, Zero)),
    gasHSource_(false),
    qrHSource_(false),
    useChemistrySolvers_(true),
    filmRegionName_("filmRegion"),
    filmRegionFound_(false),
    mCharBurnt_
    (
        IOobject
        (
            "mCharBurnt",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    charOxiQdot_
    (
        IOobject
        (
            "charOxiQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    TcharStart_(coeffs().lookupOrDefault<scalar>("TcharStart", 700)),
    TcharOxiMax_(coeffs().lookupOrDefault<scalar>("TcharOxiMax", 800)),
    charIndex_(-1),
    charring_(true)  // Will be set to false if 'char' not found in species list
{
    if (active_)
    {
        read(dict);
    }

    charIndex_ = solidThermo_->composition().species()["char"];

    if (charIndex_ < 0)
    {
        charring_ = false;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDimFF::~reactingOneDimFF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactingOneDimFF::addMassSources(const label patchi, const label facei)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchi)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][facei];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar reactingOneDimFF::solidRegionDiffNo() const
{
    scalar DiNum = -GREAT;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            sqr(regionMesh().surfaceInterpolation::deltaCoeffs())
           *fvc::interpolate(kappa())
           /fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.primitiveField())*time().deltaTValue();
    }

    return returnReduce(DiNum, maxOp<scalar>());
}


scalar reactingOneDimFF::maxDiff() const
{
    return maxDiff_;
}

label reactingOneDimFF::charIndex() const
{
    return charIndex_;
}

const volScalarField& reactingOneDimFF::rho() const
{
    return rho_;
}


const volScalarField& reactingOneDimFF::T() const
{
    return solidThermo_->T();
}


const tmp<volScalarField> reactingOneDimFF::Cp() const
{
    return solidThermo_->Cp();
}


tmp<volScalarField> reactingOneDimFF::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> reactingOneDimFF::kappa() const
{
    return solidThermo_->kappa();
}


const surfaceScalarField& reactingOneDimFF::phiGas() const
{
    return phiGas_;
}


void reactingOneDimFF::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    if(!moveMesh_)
    {
        forAll(h_, celli)
        {
            solidChemistry_->setCellReacting(celli, true);
        }
    }
    else
    {
        surfaceScalarField::Boundary& phiGasBf =
            phiGas_.boundaryFieldRef();
        label totalFaceId = 0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchi = intCoupledPatchIDs_[i];
            const scalarField& cellVol = regionMesh().V();
            scalarField& phiGasp = phiGasBf[patchi];
            forAll(phiGasp, facei)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId++];
                scalar areaF = mag(regionMesh().Sf().boundaryField()[patchi][facei]);
                forAll(cells, ci)
                {
                    scalar cellThickness = cellVol[cells[ci]]/areaF;
                    if(cellThickness <= minimumDelta_)
                    {
                        solidChemistry_->setCellReacting(cells[ci], false);
                    }
                    else
                    {
                        solidChemistry_->setCellReacting(cells[ci], true);
                    }
                }
            }
        }
    }
}


void reactingOneDimFF::postEvolveRegion()
{
    pyrolysisModel::postEvolveRegion();
}


void reactingOneDimFF::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    if (useChemistrySolvers_)
    {
        solidChemistry_->solve(time().deltaTValue());
    }
    else
    {
        solidChemistry_->calculate();
    }

    solveContinuity();

    chemistryQdot_ = solidChemistry_->Qdot()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    calculateMassTransfer();

    solidThermo_->correct();

    // Updating emissivity and absorptivity at boundaries
    updateBndEmmAbs();

    // Updating film transmissivity at boundaries
    if (filmRegionFound_)
    {
        updateFilmTransmissivity();
    }

    //luwi--enforce Temperature limits
    scalarField& T = solidThermo_->T();
    T = max(Tmin_, min(T, Tmax_));

    Info<< "pyrolysis min/max(T) = "
        << gMin(solidThermo_->T().primitiveField())
        << ", "
        << gMax(solidThermo_->T().primitiveField())
        << endl;
}


void reactingOneDimFF::info()
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info << "lib/regionModels/pyrolysisModels/reactingOneDimFF/reactingOneDimFF.C" << nl;
    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //
