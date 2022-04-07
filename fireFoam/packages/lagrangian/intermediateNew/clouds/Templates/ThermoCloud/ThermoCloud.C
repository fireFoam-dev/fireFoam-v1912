/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "ThermoCloud.H"
#include "ThermoParcel.H"

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setModels()
{
    heatTransferModel_.reset
    (
        HeatTransferModel<ThermoCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    TIntegrator_.reset
    (
        integrationScheme::New
        (
            "T",
            this->solution().integrationSchemes()
        ).ptr()
    );

    this->subModelProperties().lookup("radiation") >> radiation_;
    
//<fmglobal>
    // Confirm settings in radiationProperties file
    const IOdictionary fileRadProps
	(
	    IOobject
	    (
		"radiationProperties",
		this->mesh().time().constant(),
		this->mesh(), 
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
        );

    if (radiation_)
    {    
	//...cloudAbsorptionEmission parameters would be in model2 of binaryAbsorptionEmission dictionary
	if (
	    fileRadProps.lookupOrDefault("radiation",false)
	    && fileRadProps.lookupOrDefault<word>("radiationModel","none") == "fvDOM"
	    && fileRadProps.lookupOrDefault<word>("absorptionEmissionModel","none") == "binaryAbsorptionEmission"
	   )
	{
	    const dictionary& cloudAbsEmDict
		(
		    fileRadProps.subDict("binaryAbsorptionEmissionCoeffs").subDict("model2")
		);

	    radiation_ =
		(
		    cloudAbsEmDict.lookupOrDefault<word>("absorptionEmissionModel","none") == "cloudAbsorptionEmission"
		);
	}
    }

    if (radiation_)
    {
	const dictionary& cloudAbsEmCoeffs
	    (
		fileRadProps.subDict("binaryAbsorptionEmissionCoeffs")
		.subDict("model2")
		.subDict("cloudAbsorptionEmissionCoeffs")
	    );
	coupledRadiation_ = cloudAbsEmCoeffs.lookupOrDefault("coupledRadiation",true);
	radProp_ = cloudAbsEmCoeffs.lookupOrDefault<word>("radiationProperty","constRad");
	Switch grayModel = cloudAbsEmCoeffs.lookupOrDefault("grayModel",false);
	nBands_ = 1;
	const word scatModel(fileRadProps.lookup("scatterModel"));
	setCScat_ =
	    (
		scatModel == "cloudScatter"
             && fileRadProps.lookupOrDefault<word>("phaseFunction","preCorrected") == "preCorrected"
	    );

	// Read droplet property parameters from file
	const IOdictionary propDict
	(
	    IOobject
	    (
		"multimediaRadProperties",
		this->mesh().time().constant(),
		this->mesh(), 
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	);
	const dictionary& radProp( propDict.subDict("dropletProperties") );

        if (radProp_ == "constRad")
        {
            constAbsEff_  = radProp.lookupOrDefault<scalar>("absEfficiency",0.5);
            constSctEff_  = radProp.lookupOrDefault<scalar>("sctEfficiency",0.5);
        }
        else if (radProp_ == "diaBanded")
        {
	    List<Vector2D<scalar>> iBands;
	    (propDict.subDict("allBands")).lookup("bandLimits") >> iBands;
            beamLen_ = radProp.lookupOrDefault<scalar>("beamLength",0.2);

            radProp.lookup("energyFraction") >> energyFrac_;
	    if ((numPropBands_ = energyFrac_.size()) != iBands.size())
	    {
		FatalErrorInFunction
		    << " in file constant/multimediaRadProperties" <<nl
		    << "  Droplet radiation properties specified for "<<numPropBands_
		    << " instead of " << iBands.size()
		    << " bands\n"
		    << exit(FatalError);
	    }
            radProp.lookup("diaList") >> diaVal_;
	    label nDia = diaVal_.size();
	    absEff_.setSize(nDia, energyFrac_*0);
	    sctEff_.setSize(nDia, energyFrac_*0);
	    asyFac_.setSize(nDia, energyFrac_*0);

            radProp.lookup("absEfficiency") >> absEff_;
            radProp.lookup("sctEfficiency") >> sctEff_;
            radProp.lookup("asymmetryFactor") >> asyFac_;

	    if (!grayModel)
		nBands_ = numPropBands_;
        }
        else
        {
	    FatalErrorInFunction
		<< "Invalid keyword: "<<radProp_
		<< nl
		<< "  Available options are:\n"
		<< "  (\n"
		<< "    constRad\n"
		<< "    diaBanded\n"
		<< "  )\n"
		<< exit(FatalError);
        }

        radAreaP_.setSize(nBands_);
        radAreaPSc_.setSize(nBands_);
        radT4_.setSize(nBands_);
        radAreaPT4_.setSize(nBands_);
        radAreaPScAsy_.setSize(nBands_);
        CScat_.setSize(nBands_);
 
        scalar tInit=0.;
        CScatCoeffs_.setSize(4,tInit);
 

        forAll(radAreaP_, bandI)
        {
            radAreaP_.set
            (
                bandI,
                new volScalarField::Internal
                 (
                     IOobject
                     (
                         this->name() + ":radAreaP" + "_" + name(bandI),
                         this->db().time().timeName(),
                         this->db(),
                         IOobject::READ_IF_PRESENT,
                         IOobject::AUTO_WRITE
                     ),
                     this->mesh(),
                    dimensionedScalar(dimArea, Zero)
                 )           
            );
        }

        forAll(radAreaPSc_, bandI)
        {
            radAreaPSc_.set
            (
             bandI,
             new volScalarField::Internal
                (
                    IOobject
                    (
                        this->name() + ":radAreaPSc" + "_" + name(bandI),
                        this->db().time().timeName(),
                        this->db(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimArea, Zero)
                )
            );
        }

        forAll(radT4_, bandI)
        {
            radT4_.set
                (
                    bandI,
                    new volScalarField::Internal
                    (
                        IOobject
                        (
                            this->name() + ":radT4" + "_" + name(bandI),
                            this->db().time().timeName(),
                            this->db(),
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        this->mesh(),
                        dimensionedScalar(dimArea, Zero)
                    )
                );
        }

        forAll(radAreaPT4_, bandI)
        {
            radAreaPT4_.set
                (
                 bandI,
                 new volScalarField::Internal
                 (
                  IOobject
                  (
                   this->name() + ":radAreaPT4" + "_" + name(bandI),
                   this->db().time().timeName(),
                   this->db(),
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
                  ),
                  this->mesh(),
                  dimensionedScalar(dimArea, Zero)
                 )
                );
        }

        forAll(radAreaPScAsy_, bandI)
        {
            radAreaPScAsy_.set
                (
                 bandI,
                 new volScalarField::Internal
                 (
                  IOobject
                  (
                   this->name() + ":radAreaPScAsy" + "_" + name(bandI),
                   this->db().time().timeName(),
                   this->db(),
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
                  ),
                  this->mesh(),
                  dimensionedScalar(dimArea, Zero)
                 )
                );
        }

        forAll(CScat_, bandI)
        {
            CScat_.set
                (
                 bandI,
                 new volScalarField::Internal
                 (
                  IOobject
                  (
                   this->name() + ":CScat" + "_" + name(bandI),
                   this->db().time().timeName(),
                   this->db(),
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
                  ),
                  this->mesh(),
                  dimensionedScalar(dimArea, Zero)
                 )
                );
        }

    }
//</fmglobal>    
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::cloudReset(ThermoCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    heatTransferModel_.reset(c.heatTransferModel_.ptr());
    TIntegrator_.reset(c.TIntegrator_.ptr());

    radiation_ = c.radiation_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
    coupledRadiation_(false), // ankur
    //radModel_(T_),
    nBands_(0), // ankur
    radProp_(), // ankur
    constAbsEff_(), // ankur
    constSctEff_(), // ankur
    numPropBands_(), // ankur
    beamLen_(0.2), // ankur
    energyFrac_(), // ankur
    diaVal_(), // ankur
    absEff_(), // ankur
    sctEff_(), // ankur
    asyFac_(), // ankur
//    radAreaP_(nullptr), // ankur
//    radT4_(nullptr), // ankur
//    radAreaPT4_(nullptr), // ankur
    CScatCoeffs_(), // ankur
    setCScat_(true), // ankur
    setCScatCoeffs_(true), // ankur
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy, Zero)
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTemperature, Zero)
        )
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    ThermoCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(c.heatTransferModel_->clone()),
    TIntegrator_(c.TIntegrator_->clone()),
    radiation_(c.radiation_),
    coupledRadiation_(c.coupledRadiation_), // ankur
    //radModel_(c.radModel_), // ankur
    nBands_(c.nBands_), // ankur
    //nBands_(0), // ankur
    radProp_(c.radProp_), // ankur
    constAbsEff_(c.constAbsEff_), // ankur
    constSctEff_(c.constSctEff_), // ankur
    numPropBands_(c.numPropBands_), // ankur
    beamLen_(c.beamLen_), // ankur
    energyFrac_(c.energyFrac_), // ankur
    diaVal_(c.diaVal()), // ankur
    absEff_(c.absEff_), // ankur
    sctEff_(c.sctEff_), // ankur
    asyFac_(c.asyFac_), // ankur
//    radAreaP_(nullptr), // ankur
//    radT4_(nullptr), // ankur
//    radAreaPT4_(nullptr), // ankur
    CScatCoeffs_(c.CScatCoeffs_), // ankur
    setCScat_(c.setCScat_), // ankur
    setCScatCoeffs_(c.setCScatCoeffs_), // ankur
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsTrans()
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsCoeff()
        )
    )
{
    if (radiation_)
    {
        // ankur
        nBands_=c.nBands();
        radAreaP_(nBands_);
        radAreaPSc_(nBands_);
        radT4_(nBands_);
        radAreaPT4_(nBands_);
        radAreaPScAsy_(nBands_);
        CScat_(nBands_);

        forAll(radAreaP_, bandI)
        {
         radAreaP_.set
         (
           bandI,
           new volScalarField::Internal
            (
                IOobject
                (
                    //this->name() + ":radAreaP" + "_" + name(bandI), // kvm, name(bandI) won't compile
                    this->name() + ":radAreaP",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radAreaP(bandI)
            )           
         );
        }

        forAll(radAreaPSc_, bandI)
        {
         radAreaPSc_.set
         (
           bandI,
           new volScalarField::Internal
            (
                IOobject
                (
                    //this->name() + ":radAreaPSc" + "_" + name(bandI),
                    this->name() + ":radAreaPSc" + "_" ,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radAreaPSc(bandI)
            )
         );
        }

        forAll(radT4_, bandI)
        {
         radT4_.set
         (
           bandI,
           new volScalarField::Internal
            (
                IOobject
                (
                    //this->name() + ":radT4" + "_" + name(bandI),
                    this->name() + ":radT4" + "_",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radT4(bandI)
            )
         );
        }

        forAll(radAreaPT4_, bandI)
        {
         radAreaPT4_.set
         (
           bandI,
           new volScalarField::Internal
            (
                IOobject
                (
                    //this->name() + ":radAreaPT4" + "_" + name(bandI),
                    this->name() + ":radAreaPT4" + "_",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radAreaPT4(bandI)
            )
        );
        }

        forAll(radAreaPScAsy_, bandI)
        {
         radAreaPScAsy_.set
         (
           bandI,
           new volScalarField::Internal
            (
                IOobject
                (
                    //this->name() + ":radAreaPScAsy" + "_" + name(bandI),
                    this->name() + ":radAreaPScAsy",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.radAreaPScAsy(bandI)
            )
         );
        }

        forAll(CScat_, bandI)
        {
         CScat_.set
         (
           bandI,
           new volScalarField::Internal
            (
                IOobject
                (
                    //this->name() + ":CScat" + "_" + name(bandI),
                    this->name() + ":CScat" + "_",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.CScat(bandI)
            )
         );
        }

//        radAreaP_.set
//        (
//            new volScalarField::Internal
//            (
//                IOobject
//                (
//                    this->name() + ":radAreaP",
//                    this->db().time().timeName(),
//                    this->db(),
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE,
//                    false
//                ),
//                c.radAreaP()
//            )
//        );
//
//        radT4_.set
//        (
//            new volScalarField::Internal
//            (
//                IOobject
//                (
//                    this->name() + ":radT4",
//                    this->db().time().timeName(),
//                    this->db(),
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE,
//                    false
//                ),
//                c.radT4()
//            )
//        );
//
//        radAreaPT4_.set
//        (
//            new volScalarField::Internal
//            (
//                IOobject
//                (
//                    this->name() + ":radAreaPT4",
//                    this->db().time().timeName(),
//                    this->db(),
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE,
//                    false
//                ),
//                c.radAreaPT4()
//            )
//        );
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const ThermoCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    thermo_(c.thermo()),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
//    radAreaP_(nullptr),
//    radT4_(nullptr),
//    radAreaPT4_(nullptr),
    nBands_(0),
    hsTrans_(nullptr),
    hsCoeff_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.T() = constProps_.T0();
    parcel.Cp() = constProps_.Cp0();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ThermoCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    hsTrans_->field() = 0.0;
    hsCoeff_->field() = 0.0;

    if (radiation_)
    {
        // ankur
        forAll(radAreaP_, bandI)
        {
            radAreaP_[bandI].field() = 0.0;
            radAreaPSc_[bandI].field() = 0.0;
            radT4_[bandI].field() = 0.0;
            radAreaPT4_[bandI].field() = 0.0;
            radAreaPScAsy_[bandI].field() = 0.0;
            CScat_[bandI].field() = 0.0;
        }
//        setCScat_ = false;
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::relaxSources
(
    const ThermoCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(hsTrans_(), cloudOldTime.hsTrans(), "h");
    this->relax(hsCoeff_(), cloudOldTime.hsCoeff(), "h");

    if (radiation_)
    {
        // ankur
        forAll(radAreaP_, bandI)
        {
            this->relax(radAreaP_[bandI], cloudOldTime.radAreaP(bandI), "radiation");
            this->relax(radAreaPSc_[bandI], cloudOldTime.radAreaPSc(bandI), "radiation");
            this->relax(radT4_[bandI], cloudOldTime.radT4(bandI), "radiation");
            this->relax(radAreaPT4_[bandI], cloudOldTime.radAreaPT4(bandI), "radiation");
            this->relax(radAreaPScAsy_[bandI], cloudOldTime.radAreaPScAsy(bandI), "radiation");
            this->relax(CScat_[bandI], cloudOldTime.CScat(bandI), "radiation");
        }
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(hsTrans_(), "h");
    this->scale(hsCoeff_(), "h");

    if (radiation_)
    {
        // ankur
        forAll(radAreaP_, bandI)
        {
            this->scale(radAreaP_[bandI], "radiation");
            this->scale(radAreaPSc_[bandI], "radiation");
            this->scale(radT4_[bandI], "radiation");
            this->scale(radAreaPT4_[bandI], "radiation");
            this->scale(radAreaPScAsy_[bandI], "radiation");
            this->scale(CScat_[bandI], "radiation");
        }
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::preEvolve()
{
    CloudType::preEvolve();

    this->pAmbient() = thermo_.thermo().p().average().value();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}

template<class CloudType>
void Foam::ThermoCloud<CloudType>::updateProp() // ankur
{

    if (this->solution().active() && (radiation_) && setCScat_ )
    {
        if (setCScatCoeffs_) 
        {
            setCScatCoeffs();
            setCScatCoeffs_ = false;
        }

        setCScat();
    }

}

template<class CloudType>
void Foam::ThermoCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::info()
{
    CloudType::info();

    Info<< "    Temperature min/max             = " << Tmin() << ", " << Tmax()
        << endl;
}


// ************************************************************************* //
