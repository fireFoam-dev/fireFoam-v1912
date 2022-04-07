/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "fvDOM.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "unitConversion.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM, 0);
        addToRadiationRunTimeSelectionTables(fvDOM);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOM::rotateInitialRays(const vector& sunDir)
{
    // Rotate Y spherical cordinates to Sun direction.
    // Solid angles on the equator are better fit for planar radiation
    const tensor coordRot = rotationTensor(vector(0, 1, 0), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }
}


void Foam::radiation::fvDOM:: alignClosestRayToSun(const vector& sunDir)
{
    label SunRayId(-1);
    scalar maxSunRay = -GREAT;

    // Looking for the ray closest to the Sun direction
    forAll(IRay_, rayId)
    {
        const vector& iD = IRay_[rayId].d();
        scalar dir = sunDir & iD;
        if (dir > maxSunRay)
        {
            maxSunRay = dir;
            SunRayId = rayId;
        }
    }

    // Second rotation to align colimated radiation with the closest ray
    const tensor coordRot = rotationTensor(IRay_[SunRayId].d(), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }

    Info << "Sun direction : " << sunDir << nl << endl;
    Info << "Sun ray ID : " << SunRayId << nl << endl;
}


void Foam::radiation::fvDOM::updateRaysDir()
{
    solarCalculator_->correctSunDirection();
    const vector sunDir = solarCalculator_->direction();

    // First iteration
    if (updateTimeIndex_ == 0)
    {
        rotateInitialRays(sunDir);
        alignClosestRayToSun(sunDir);
    }
    else if (updateTimeIndex_ > 0)
    {
        alignClosestRayToSun(sunDir);
    }
}


void Foam::radiation::fvDOM::initialise()
{
    coeffs_.readIfPresent("useExternalBeam", useExternalBeam_);

    if (useExternalBeam_)
    {
        coeffs_.readEntry("spectralDistribution", spectralDistribution_);

        spectralDistribution_ =
            spectralDistribution_/sum(spectralDistribution_);

        const dictionary& solarDict = this->subDict("solarCalculatorCoeffs");
        solarCalculator_.reset(new solarCalculator(solarDict, mesh_));

        if (mesh_.nSolutionD() != 3)
        {
            FatalErrorInFunction
                << "External beam model only available in 3D meshes "
                << abort(FatalError);
        }

        if (solarCalculator_->diffuseSolarRad() > 0)
        {
            FatalErrorInFunction
                << "External beam model does not support Diffuse "
                << "Solar Radiation. Set diffuseSolarRad to zero"
                << abort(FatalError);
        }
        if (spectralDistribution_.size() != nLambda_)
        {
            FatalErrorInFunction
                << "The epectral energy distribution has different bands "
                << "than the absoprtivity model "
                << abort(FatalError);
        }
    }

    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        nRay_ = 4*nPhi_*nTheta_;

        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi/(2.0*nPhi_);
        scalar deltaTheta = pi/nTheta_; // ankur
        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0; // ankur
            // ankur
            if (uniSolidAngles_)
            {
              scalar muBeg = 1. - (n - 1.)*2./nTheta_;   // mu is cos(theta)
              scalar muEnd = 1. - n*(2./nTheta_);
              scalar thetaBeg = Foam::acos(muBeg);
              scalar thetaEnd = Foam::acos(muEnd);
              deltaTheta = thetaEnd - thetaBeg;
              thetai = thetaBeg + (deltaTheta/2.); 
            }
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        *absorptionEmission_,
                        scatter_, // ankur
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 4*nPhi_;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi/(2.0*nPhi_);
        label i = 0;
        for (label m = 1; m <= 4*nPhi_; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    scatter_, // ankur
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }
    // 1D
    else
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 2;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    scatter_, // ankur
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    // ankur, construct incident radiation for each band
    forAll(GLambda_, lambdaI)
    {
        GLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "GLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                G_
            )
        );
    }

    // ankur, construct energy fraction for each band
    forAll(enFracLambda_, lambdaI)
    {
        enFracLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "enFracLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }
    
    // luwi
    // Construct boundary fluxes for each band
    // --incident radiation
    forAll(qinLambda_, lambdaI)
    {	
        qinLambda_.set
    	(
	        lambdaI,
	        new volScalarField
	        (
	        	IOobject
	        	(
	        	    "qinLambda_" + Foam::name(lambdaI) ,
	        	    mesh_.time().timeName(),
	        	    mesh_,
	        	    IOobject::NO_READ,
	        	    IOobject::NO_WRITE
	        	),
        		qin_
	        )
    	);
    }
    // luwi
    // --emitted radiation
    forAll(qemLambda_, lambdaI)
    {
    	qemLambda_.set
    	(
	        lambdaI,
	        new volScalarField
	        (
	            IOobject
	        	(
	                "qemLambda_" + Foam::name(lambdaI) ,
	                mesh_.time().timeName(),
	        	    mesh_,
	        	    IOobject::NO_READ,
	        	    IOobject::NO_WRITE
	        	),
        		qem_
	        )
	    );
    }
    // luwi
    // --net radiation
    forAll(qrLambda_, lambdaI)
    {
    	qrLambda_.set
    	(
	        lambdaI,
	        new volScalarField
	        (
	            IOobject
	        	(
	                "qrLambda_" + Foam::name(lambdaI) ,
	                mesh_.time().timeName(),
	        	    mesh_,
	        	    IOobject::NO_READ,
	        	    IOobject::NO_WRITE
	        	),
        		qr_
	        )
	    );
    }

    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    if (useExternalBeam_)
    {
        // Rotate rays for Sun direction
        updateRaysDir();
    }

    scalar totalOmega = 0;
    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        totalOmega += IRay_[rayId].omega();
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "dAve : "
            << '\t' << IRay_[rayId].dAve() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << " : " << "d : "
            << '\t' << IRay_[rayId].d() << nl;
    }

    Info << "Total omega : " << totalOmega << endl;

    Info<< endl;

//<fmglobal>
    if (this->found("emissionTemperature") )
    {
	const dictionary& emissionTDict = this->subDict("emissionTemperature");
	Ctri_  = emissionTDict.lookupOrDefault<scalar>("C_tri",1.25);
	Ctvar_ = emissionTDict.lookupOrDefault<scalar>("C_tvar",2.0);
	TRI_ = emissionTDict.lookupOrDefault("TRI",false);
	if (TRI_)
	{
	    forAll(delta2_, i)
	    {	
	        delta2_[i] = pow(mesh_.V()[i],(2./3.));
	    }
	}
    }
    else
    {
        TRI_ = false;
    }
//</fmglobal>

    coeffs_.readIfPresent("useSolarLoad", useSolarLoad_);

    if (useSolarLoad_)
    {
        if (useExternalBeam_)
        {
            FatalErrorInFunction
                << "External beam with fvDOM can not be used "
                << "with the solar load model"
                << abort(FatalError);
        }
        const dictionary& solarDict = this->subDict("solarLoadCoeffs");
        solarLoad_.reset(new solarLoad(solarDict, T_));

        if (solarLoad_->nBands() != this->nBands())
        {
            FatalErrorInFunction
                << "Requested solar radiation with fvDOM. Using "
                << "different number of bands for the solar load is not allowed"
                << abort(FatalError);
        }

        Info<< "Creating Solar Load Model " << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::fvDOM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    GLambda_(nLambda_), // ankur
    enFracLambda_(nLambda_), // ankur
    qrLambda_(nLambda_), // luwi
    qinLambda_(nLambda_), // luwi
    qemLambda_(nLambda_), // luwi
//<fmglobal>
    delta2_
    (
	IOobject
	(
	    "delta2",
	    mesh_.time().timeName(),
	    mesh_,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh_,
	dimensionedScalar(dimLength*dimLength, Zero)
    ),
//</fmglobal>
    blackBody_(nLambda_, T),
    IRay_(0),
    tolerance_
    (
        coeffs_.lookupOrDefaultCompat("tolerance", {{"convergence", 1712}}, 0.0)
    ),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    uniSolidAngles_(coeffs_.lookupOrDefault<Switch>("uniformSolidAngles",false)), // ankur
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)
    ),
    useExternalBeam_(false),
    spectralDistribution_(),
    solarCalculator_(),
    updateTimeIndex_(0)
{
    initialise();
}


Foam::radiation::fvDOM::fvDOM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    GLambda_(nLambda_), // ankur
    enFracLambda_(nLambda_), // ankur
    qrLambda_(nLambda_), // luwi
    qinLambda_(nLambda_), // luwi
    qemLambda_(nLambda_), // luwi
//<fmglobal>
    delta2_
    (
	IOobject
	(
	    "delta2",
	    mesh_.time().timeName(),
	    mesh_,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh_,
	dimensionedScalar(dimLength*dimLength, Zero)
    ),
//</fmglobal>
    blackBody_(nLambda_, T),
    IRay_(0),
    tolerance_
    (
        coeffs_.lookupOrDefaultCompat("tolerance", {{"convergence", 1712}}, 0.0)
    ),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    uniSolidAngles_(coeffs_.lookupOrDefault<Switch>("uniformSolidAngles",false)), // ankur
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)
    ),
    useExternalBeam_(false),
    spectralDistribution_(),
    solarCalculator_(),
    updateTimeIndex_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresentCompat
        (
            "tolerance", {{"convergence", 1712}}, tolerance_
        );
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }

    return false;
}


void Foam::radiation::fvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    absorptionEmission_->correctEnFrac(enFracLambda_, blackBody_); // ankur, to update energy fraction 
    
    //- Efficient computation of normalized Henyey-Greenstein phase functions...
    //  compute and store all unique values once
    scatter_->updatePhaseFuncs(*this); // luwi

    // updateBlackBodyEmission();  // ankur, this not needed now, since the above call updates the fraction.. Also, the data for storing spectral fraction has now been defined in fvDOM class itself now..
    
    if (useSolarLoad_)
    {
        solarLoad_->calculate();
    }

    if (useExternalBeam_)
    {
        switch (solarCalculator_->sunDirectionModel())
        {
            case solarCalculator::mSunDirConstant:
            {
                break;
            }
            case solarCalculator::mSunDirTracking:
            {
                label updateIndex = label
                (
                    mesh_.time().value()
                   /solarCalculator_->sunTrackingUpdateInterval()
                );

                if (updateIndex > updateTimeIndex_)
                {
                    Info << "Updating Sun position..." << endl;
                    updateTimeIndex_ = updateIndex;
                    updateRaysDir();
                }
                break;
            }
        }
    }

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0;
    label radIter = 0;
    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;

        radIter++;
        maxResidual = 0;
        forAll(IRay_, rayI)
        {
//<fmglobal>		
	    if (radIter==1)
	    {
		IRay_[rayI].resetConvLambda();
	    }
	    //luwi note: in-scattering couples all rays; revisit ray convergence criterion below
//</fmglobal>
            if (!rayIdConv[rayI])
            {
                scalar maxBandResidual = IRay_[rayI].correct();
                maxResidual = max(maxBandResidual, maxResidual);

                if (maxBandResidual < tolerance_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }

    } while (maxResidual > tolerance_ && radIter < maxIter_);

    updateG();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            // Only include continuous phase emission
            //4*absorptionEmission_->aCont()*physicoChemical::sigma  // ankur, commenting this line.. not correct for wide-band/WSGG models
            4*this->aCont()*physicoChemical::sigma // luwi
        )
    );
}


Foam::tmp<Foam::volScalarField::Internal> // ankur
Foam::radiation::fvDOM::Ru() const
{

//    const volScalarField::Internal& G =
//        G_();
//
//    const volScalarField::Internal E =
//        absorptionEmission_->ECont()()();
//
//    // Only include continuous phase absorption
//    const volScalarField::Internal a =
//        absorptionEmission_->aCont()()();
//
//    return a*G - E;

    // ankur
    tmp<volScalarField::Internal> tRu
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tRu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), Zero)
        )
    );

    // ankur
    for (label j=0; j < nLambda_; j++)
    {
        //Info << "j = " << j << nl;
        //Info << "aCont = " << absorptionEmission_->aCont(j)()() << nl;
        //Info << "GLambda = " << GLambda_[j]() << nl;
        //Info << "Econt = " << absorptionEmission_->ECont(j)()() << nl;
        tRu.ref() += absorptionEmission_->aCont(j)()()*GLambda_[j]() - absorptionEmission_->ECont(j)()();
    }

    return tRu; // ankur

}


void Foam::radiation::fvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}

// luwi
// continuous phase absorption coefficient
Foam::tmp<Foam::volScalarField>
Foam::radiation::fvDOM::aCont() const
{
    tmp<volScalarField> taCont
    (
        new volScalarField
        (
            IOobject
            (
                "taCont",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );

    for (label j=0; j < nLambda_; j++)
    {
        taCont.ref() += absorptionEmission_->aCont(j)*enFracLambda_[j]; 
    }   

    return taCont;
}


void Foam::radiation::fvDOM::updateG()
{
    G_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qr_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qem_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qin_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);

    // luwi
    if ( nLambda_ > 1 )
    {
        // non-grey
        forAll(GLambda_,iLambda)
        {
        	GLambda_[iLambda] = dimensionedScalar(dimMass/pow3(dimTime), Zero);
        	qinLambda_[iLambda] = dimensionedScalar(dimMass/pow3(dimTime), Zero);
        	qemLambda_[iLambda] = dimensionedScalar(dimMass/pow3(dimTime), Zero);
        	qrLambda_[iLambda]  = dimensionedScalar(dimMass/pow3(dimTime), Zero);
        }

        forAll(IRay_, rayI)
        {
            IRay_[rayI].addIntensity();
            G_ += IRay_[rayI].I()*IRay_[rayI].omega();
            
            forAll(GLambda_,iLambda)
            {
                GLambda_[iLambda] += IRay_[rayI].ILambda(iLambda)*IRay_[rayI].omega();
                // update boundary spectral heat fluxes
                qinLambda_[iLambda].boundaryFieldRef() += IRay_[rayI].qinLambda(iLambda).boundaryField();
                qemLambda_[iLambda].boundaryFieldRef() += IRay_[rayI].qemLambda(iLambda).boundaryField();
                qrLambda_[iLambda].boundaryFieldRef()  += IRay_[rayI].qrLambda(iLambda).boundaryField();
            }
        }
        // update boundary total heat fluxes
        forAll(GLambda_,iLambda)
        {
            qin_.boundaryFieldRef() += qinLambda_[iLambda].boundaryFieldRef();
            qem_.boundaryFieldRef() += qemLambda_[iLambda].boundaryFieldRef();
            qr_.boundaryFieldRef()  += qrLambda_[iLambda].boundaryFieldRef();
        }
    }
    else
    {
        // grey
        forAll(IRay_, rayI)
        {            
	    IRay_[rayI].addIntensity();
            G_ += IRay_[rayI].I()*IRay_[rayI].omega();
            GLambda_[0] += IRay_[rayI].ILambda(0)*IRay_[rayI].omega();      // pc
   //Info << "IRay_ " << IRay_[rayI].I() <<  nl; 
   //Info << "omega " << IRay_[rayI].omega() <<  nl; 
            // update boundary heat fluxes
            qin_.boundaryFieldRef() += IRay_[rayI].qin().boundaryField();
            qem_.boundaryFieldRef() += IRay_[rayI].qem().boundaryField();
            qr_.boundaryFieldRef()  += IRay_[rayI].qr().boundaryField();
        }
	GLambda_[0] = G_;
    }
}


void Foam::radiation::fvDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // Assuming name is in the form: CHARS_rayId_lambdaId
    const auto i1 = name.find('_');
    const auto i2 = name.find('_', i1+1);

    rayId    = readLabel(name.substr(i1+1, i2-i1-1));
    lambdaId = readLabel(name.substr(i2+1));
}

// ankur
Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::inScatEnergy(const label iLambda, const label sourDir) const
{
    tmp<volScalarField> tInScatEn
    (
        new volScalarField
        (
            IOobject
            (
                "inScatEn",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );

    forAll(IRay_, rayI)
    {
        tInScatEn.ref() += scatter_->pFunc(iLambda,rayI,sourDir)*IRay_[rayI].ILambda(iLambda)*IRay_[rayI].omega();

        //if (rayI != sourDir)
        //{
           //dimensionedScalar temp1(dimLength, 1.0);
           //tInScatEn.ref() += scatter_->sigmaEff(iLambda)*IRay_[rayI].ILambda(iLambda)*IRay_[rayI].omega()/(4.*pi - IRay_[rayI].omega());
        //}
//Info << "scatIn: " << scatter_->pFunc(iLambda,rayI,sourDir).ref().primitiveFieldRef() << endl;
    }

    if( !scatter_->cloudCorrectedSigma() )
	tInScatEn.ref() *=
	    (
            dimensionedScalar("one",dimLength,1.0)
            *scatter_->sigmaEff(iLambda,-1,-1.)/(4.*pi)
	    );  

    return tInScatEn;
}
//<fmglobal>
Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::updateT4fac()
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

    if( TRI_ )
    {
	T4fac_ = 1. + 6.*Ctri_*Ctvar_*delta2_*pow(mag(fvc::grad(T_)),2)/(T_*T_);    
	tT4fac.ref() = T4fac_;
    }

    return tT4fac;
}
//</fmglobal>

const Foam::solarCalculator& Foam::radiation::fvDOM::solarCalc() const
{
    return solarCalculator_();
}


// ************************************************************************* //
