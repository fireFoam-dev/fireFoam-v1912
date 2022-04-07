/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "boxFixedBandsAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(boxFixedBandsAbsorptionEmission, 0);
        addToRunTimeSelectionTable(absorptionEmissionModel, boxFixedBandsAbsorptionEmission, dictionary);
    }
}

namespace Foam
{
    const Foam::Enum
    <
        Foam::radiation::boxFixedBandsAbsorptionEmission::absCoeffType
    >abscoeffTypeNames_
    {
        { Foam::radiation::boxFixedBandsAbsorptionEmission::absCoeffType::absCoeffAverage, "average" },
        { Foam::radiation::boxFixedBandsAbsorptionEmission::absCoeffType::absCoeffLocal, "local" },
        { Foam::radiation::boxFixedBandsAbsorptionEmission::absCoeffType::absCoeffAverageAlpha, "averageAlpha"}
    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::boxFixedBandsAbsorptionEmission::boxFixedBandsAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    boxFileName_( coeffsDict_.lookupOrDefault<word>("lookupFile","multimediaRadProperties") ),
    boxSpeciesDict_( boxDict(mesh,"boxRadiatingSpecies") ),
    boxDict_( boxDict(mesh,"boxRadiatingBands") ),
    gasThermo_(	dynamic_cast<const reactingMixture<gasHThermoPhysics>&>	( thermoRef(mesh) ) ),
    p_( thermoRef(mesh).p() ),
    T_( thermoRef(mesh).T() ),
    Y_( thermoRef(mesh).composition().Y() ),
    mixNspecies_( Y_.size() ),
    indexCO2_( gasThermo_.species()["CO2"] ),
    indexH2O_( gasThermo_.species()["H2O"] ),
    molWtCO2_( gasThermo_.speciesData()[indexCO2_].W() ),
    molWtH2O_( gasThermo_.speciesData()[indexH2O_].W() ),
    Csoot_(coeffsDict_.lookupOrDefault<scalar>("Csoot",0.0)),
    mixtureFractionSoot_( Csoot_==0? false : true ),
    T_threshold_( coeffsDict_.lookupOrDefault<scalar>("T_threshold", 1.0) ),
    method_(abscoeffTypeNames_[coeffsDict_.lookupOrDefault<word>("abscoeffType","average")])

{
    // Read list of all band limits
    (boxDict(mesh,"allBands")).lookup("bandLimits") >> iBands_;

    // Read radiating species and their spectrally-invariant properties
    radNspecies_=0;
    forAllConstIter(dictionary, boxSpeciesDict_, iter)
    {
	if ( iter().isDict() )
	{
	    const dictionary& speciesDict = iter().dict();
	    const word& key = iter().keyword();
	    //for whatever reason, label m(dict.lookup("m")) wouldn't compile
	    label m;
	    speciesDict.lookup("m") >> m;
	    const scalarList g_k = speciesDict.lookup("g_k");

	    if ( key == "H2O" )
	    {
		if ( !(m==3 && g_k==scalarList{1,1,1}) )
		{
		    FatalErrorInFunction
			<< " H2O dictionary in "
			<< boxFileName_
			<< " must have m=3 and g_k= (1 1 1)"
			<< exit(FatalError);
		}
		speciesDict.lookup("eta_k") >> etakH2O_;
		radNspecies_++;
	    }
	    else if ( key=="CO2" )
	    {
		if ( !(m==3 && g_k==scalarList{1,2,1}) )
		{
		    FatalErrorInFunction
			<< " CO2 dictionary in "
			<< boxFileName_
			<< " must have m=3 and g_k= (1 1 1)"
			<< exit(FatalError);
		}
		speciesDict.lookup("eta_k") >> etakCO2_;
		radNspecies_++;
	    }
	    else
	    {
		FatalErrorInFunction
		    << "  Invalid radiating species "
		    << key <<" found in " << boxFileName_
		    << nl
		    << "  Available options are:\n"
		    << "  (\n"
		    << "    H2O\n"
		    << "    CO2\n"
		    << "  )\n"
		    << exit(FatalError);
	    }
	}
    }

    // Read spectral bands and properties
    nBands_ = iBands_.size();
    coeffs_.setSize(nBands_);
    transparentBand_.setSize(nBands_, label(1));
    deltaEta_.setSize(nBands_);
    forAllConstIter(dictionary, boxDict_, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }

	const label maxBandID = 50;
	const word& key = iter().keyword();

	for(label iband=0; iband<maxBandID; iband++)
	{
	    if( key == "band"+name(iband) )
	    {
		const dictionary& bandDict = iter().dict();
		Vector2D<scalar> iBand;
		bandDict.lookup("bandLimits") >> iBand;
		label nBand = -1;
		forAll(iBands_, i)
		{
		    if ( iBand == iBands_[i] )
		    {
			nBand = i;
			break;
		    }
		}
		//
		if (nBand>0)
		{
		    coeffs_[nBand].initialiseEWB(bandDict);
		    //Bandwidth [ cm^-1 ]
		    deltaEta_[nBand] = (1./iBands_[nBand][0]-1./iBands_[nBand][1])*1.e-2;
		    transparentBand_[nBand] = 0;
		    break;
		}
		else if (nBand==0)
		{
		    FatalErrorInFunction
			<< key <<" in " << boxFileName_ <<".allBands is not transparent\n"
			<<" First band in spectrum must be transparent"
			<< exit(FatalError);
		}
		else
		{
		    FatalErrorInFunction
			<< key <<" not found in " << boxFileName_ <<".allBands"
			<< nl
			<< "  Available bands are:\n"
			<< iBands_
			<< exit(FatalError);
		}
	    }
	}
    }

}

const Foam::dictionary
Foam::radiation::boxFixedBandsAbsorptionEmission::boxDict(const fvMesh& mesh, const word key)
{
    // Read box model parameters from lookup file
    IOdictionary boxFileDict
    (
	IOobject
	(
	    boxFileName_,
	    mesh.time().constant(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	)
    );

    return( boxFileDict.subDict(key) );
}

const Foam::psiReactionThermo&
Foam::radiation::boxFixedBandsAbsorptionEmission::thermoRef(const fvMesh& mesh) const
{
    return ( mesh.lookupObject<psiReactionThermo>("thermophysicalProperties") );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::boxFixedBandsAbsorptionEmission::~boxFixedBandsAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::boxFixedBandsAbsorptionEmission::aCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );
    scalarField& a = ta.ref().primitiveFieldRef();
    static scalar T_avg=0.;
    static scalar rhoCO2_avg = 0.;
    static scalar rhoH2O_avg = 0.;
    static scalarList rhoCO2(Y_[0].size(),0.);
    static scalarList rhoH2O(Y_[0].size(),0.);

    if ( bandI==0 )// Compute and store spectral-band-invariant values
    {
	// Domain averages...
	if ( method_ != absCoeffLocal )
	{
	    // create temperature list to be averaged
	    volScalarField T_hold(T_);
	    volScalarField Y_mean(Y_[0]*0.);
	    scalarList YTp_avg(mixNspecies_+2,0.);
	    scalar& Temp_avg = YTp_avg[mixNspecies_];
	    scalar& p_avg    = YTp_avg[mixNspecies_+1];
	    label nT = 0;

	    forAll(T_, i)
	    {
		if( T_[i] > T_threshold_ ) // if temperature > threshold, considered within flame zone
		{
		    Temp_avg   += T_[i];
		    nT++;
		    if ( method_==absCoeffAverage )
		    {
			p_avg      += p_[i];
			for (label k=0; k<mixNspecies_; k++)
			{
			    YTp_avg[k] += Y_[k][i];
			}
		    }
		}
	    }

	    // Global averages
	    reduce(nT, sumOp<label>());
	    gSum(YTp_avg);
	    YTp_avg = YTp_avg/double(nT);
	    T_avg = Temp_avg;

	    if ( method_==absCoeffAverage )
	    {
		scalar mixMolWt = 0.;
		for (label specieI=0; specieI<mixNspecies_; specieI++)
		{
		    const scalar& molWti = gasThermo_.speciesData()[specieI].W();
		    mixMolWt += YTp_avg[specieI]/molWti;
		}
		mixMolWt = 1./mixMolWt;
		// radiating species mass densities in g/m3
		rhoCO2_avg = YTp_avg[indexCO2_]*mixMolWt*p_avg/(gasConstant_*T_avg);
		rhoH2O_avg = YTp_avg[indexH2O_]*mixMolWt*p_avg/(gasConstant_*T_avg);
	    }
	}
        // Local species densities...
	if ( method_ != absCoeffAverage )
	{
	    volScalarField mixMolWt((Y_[0]-Y_[0]));

	    forAll(Y_,specieI)
	    {
		const scalar& molWti = gasThermo_.speciesData()[specieI].W();
		mixMolWt += Y_[specieI]/molWti;
	    }
	    mixMolWt = 1/mixMolWt;

	    rhoCO2 = Y_[indexCO2_]*mixMolWt*p_/(gasConstant_*T_);
	    rhoH2O = Y_[indexH2O_]*mixMolWt*p_/(gasConstant_*T_);
	}
    }

    // Compute kappa
    if(!transparentBand_[bandI])
    {
	// rhoAlpha         [=] cm^-1/m
	scalar alpha_i = 0;
	scalar overlapAlpha_i = 0;
	switch(method_)
	{
	case absCoeffAverage:
	{
	    alpha_i = alpha(T_avg, bandI);
	    overlapAlpha_i = alpha(T_avg, -1);
	    const scalar a_avg = rhoAlpha(rhoH2O_avg,rhoCO2_avg,alpha_i,overlapAlpha_i,bandI)/deltaEta_[bandI];
	    a = a_avg;
	}
	break;

	case absCoeffLocal:
	{
	    forAll(a, cellI)
	    {
		const scalar& T_cell = T_[cellI];
		const scalar& rhoH2O_cell = rhoH2O[cellI];
		const scalar& rhoCO2_cell = rhoCO2[cellI];

		alpha_i = alpha(T_cell, bandI);
		overlapAlpha_i = alpha(T_cell, -1);
		a[cellI] = rhoAlpha(rhoH2O_cell,rhoCO2_cell,alpha_i,overlapAlpha_i,bandI)/deltaEta_[bandI];
	    }
	}
	break;

	case absCoeffAverageAlpha:
	{
	    // average density-absorption-coefficient (alpha), with local species densities (rho)
	    alpha_i = alpha(T_avg, bandI);
	    overlapAlpha_i = alpha(T_avg, -1);
	    forAll(a, cellI)
	    {
		const scalar& rhoH2O_cell = rhoH2O[cellI];
		const scalar& rhoCO2_cell = rhoCO2[cellI];

		a[cellI] = rhoAlpha(rhoH2O_cell,rhoCO2_cell,alpha_i,overlapAlpha_i,bandI)/deltaEta_[bandI];
	    }
	}
	break;

	}
    }


#if 0
    // Soot term
    if( mixtureFractionSoot_ )
    {
	const volScalarField& soot = mesh_.lookupObject<volScalarField>("soot");
	a += enFrac[bandI]*Csoot_*soot/10*T/dimensionedScalar("unity",dimTemperature,1.0)/dimensionedScalar("unity",dimLength,1.0);
    }
#endif

    return ta;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::boxFixedBandsAbsorptionEmission::eCont(const label bandI) const
{
    return aCont(bandI);
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::boxFixedBandsAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    return E;
}


inline Foam::scalar
Foam::radiation::boxFixedBandsAbsorptionEmission::alphaFactor
(
    const label eta_c,
    const word speciesName,
    const scalar T_avg
) const
{
    //S_k = 1/(1 - exp(-u_k)); where u_k(T) = C2_*eta_k/T
    scalar psiFactor = 0.;
    scalar psiScaled = 0.;

    switch (eta_c)
    {
    case 7250:// non-overlapping H2O band at 1.38 microns;  psi = S1*S3
    {
	const scalarList& eta_k = etakH2O_;
	scalar u1     = C2_*eta_k[0]/T_avg;
	scalar u3     = C2_*eta_k[2]/T_avg;
	scalar u      = u1 + u3;
	scalar uref1  = C2_*eta_k[0]/T0_;
	scalar uref3  = C2_*eta_k[2]/T0_;
	scalar uref   = uref1 + uref3;

	psiFactor = ( 1. - exp(-u) ) / ( 1. - exp(-uref) );
	psiScaled = ( (1.-exp(-uref1))*(1.-exp(-uref3)) )
                  / ( (1.-exp(-u1))*(1.-exp(-u3)) );
    }
    break;

    case 3760:
    {// non-overlapping CO2 band at 2.7 microns; psi = S1*S3
	const scalarList& eta_k = etakCO2_;
	scalar u1    = C2_*eta_k[0]/T_avg;
	scalar u3    = C2_*eta_k[2]/T_avg;
	scalar u     = u1 + u3;
	scalar uref1 = C2_*eta_k[0]/T0_;
	scalar uref3 = C2_*eta_k[2]/T0_;
	scalar uref  = uref1 + uref3;

	psiFactor = ( 1. - exp(-u) ) / ( 1. - exp(-uref) );
	psiScaled = ( (1.-exp(-uref1))*(1.-exp(-uref3)) )
	    / ( (1.-exp(-u1))*(1.-exp(-u3)) );
    }
    break;

    case 5350:
    {
	if(speciesName=="H2O")
	{// H2O band at 1.87 microns; psi = S2*S3
	    const scalarList& eta_k = etakH2O_;
	    scalar u2    = C2_*eta_k[1]/T_avg;
	    scalar u3    = C2_*eta_k[2]/T_avg;
	    scalar u     = u2 + u3;
	    scalar uref2 = C2_*eta_k[1]/T0_;
	    scalar uref3 = C2_*eta_k[2]/T0_;
	    scalar uref  = uref2 + uref3;

	    psiFactor = ( 1. - exp(-u) ) / ( 1. - exp(-uref) );
	    psiScaled = ( (1.-exp(-uref2))*(1.-exp(-uref3)) )
                      / ( (1.-exp(-u2))*(1.-exp(-u3)) );
	}
	else if(speciesName=="CO2")
	{// CO2 band at 2.0 microns (5200 cm^-1); psi = 2*S1^2*S3
	    const scalarList& eta_k = etakCO2_;
	    scalar u1    = C2_*eta_k[0]/T_avg;
	    scalar u3    = C2_*eta_k[2]/T_avg;
	    scalar u     = 2.*u1 + u3;
	    scalar uref1 = C2_*eta_k[0]/T0_;
	    scalar uref3 = C2_*eta_k[2]/T0_;
	    scalar uref  = 2.*uref1 + uref3;

	    psiFactor = ( 1. - exp(-u) ) / ( 1. - exp(-uref) );
	    psiScaled = ( (1.-exp(-uref1))*(1.-exp(-uref1))*(1.-exp(-uref3)) )
                      / ( (1.-exp(-u1))*(1.-exp(-u1))*(1.-exp(-u3)) );
	}
    }
    break;

    default:
	WarningInFunction
	    << "  Unrecognized value of eta_c:" << eta_c
	    <<nl;
	break;
    }

    return ( psiFactor*psiScaled );
}


inline Foam::scalar
Foam::radiation::boxFixedBandsAbsorptionEmission::
alpha(const scalar Temp, const label bandI) const
{

    scalar alpha = 0;
    static scalar overlapAlpha = 0;
    if ( bandI<0 )
    {// this is just a call to fetch overlapAlpha that was previously computed
	return overlapAlpha;
    }

    // Compute
    label eta_c = int(coeffs_[bandI].eta_c()); // cm-1
    const scalar& alpha_0 = coeffs_[bandI].alpha_0(); // cm-1/(g/m2)
    const scalar& overlapAlpha_0 = coeffs_[bandI].overlapAlpha_0(); // cm-1/(g/m2)
    overlapAlpha = 0;
    switch (eta_c)
    {
    case 140:// H2O rotational band at 71 microns
	alpha = alpha_0*exp(-9.*(sqrt(T0_/Temp) - 1.));
	break;

    case 1600:// non-overlapping H2O band at 6.3 microns
	alpha = alpha_0;
	break;

    case 667: // non-overlapping CO2 band at 15  microns
    case 2410:// non-overlapping CO2 band at 4.3 microns
	alpha = alpha_0;
	break;

    case 3760:// H2O band at 2.7 microns, consisting of 3 overlapping sub-bands; with one overlapping CO2 band
    {
	alpha = alpha_0;
	if( overlapAlpha_0>0 )
	    overlapAlpha = overlapAlpha_0*alphaFactor(eta_c, "CO2", Temp);
    }
    break;

    case 7250:// non-overlapping H2O band at 1.38 microns
	alpha = alpha_0*alphaFactor(eta_c, "H2O", Temp);
	break;

    case 5350:// H2O+CO2 band; 2 overlapping bands: one H2O band at 1.87 microns and one CO2 band at 2.0(?) microns
    {
	alpha = alpha_0*alphaFactor(eta_c, "H2O", Temp);
	if( overlapAlpha_0>0 )
	    overlapAlpha = overlapAlpha_0*alphaFactor(eta_c, "CO2", Temp);
    }
    break;

    default:
	WarningInFunction
	    << "  Unrecognized value of eta_c:" << eta_c
	    <<nl;
	break;
    }

    return alpha;
}


inline Foam::scalar
Foam::radiation::boxFixedBandsAbsorptionEmission::
rhoAlpha(const scalar rhoH2O, const scalar rhoCO2, const scalar alpha,
	 const scalar overlapAlpha, const label bandI) const
{
    label eta_c = int(coeffs_[bandI].eta_c()); // cm-1
    scalar rhoalpha = 0;

    switch (eta_c)
    {
    case 140: // H2O rotational band at 71 microns
    case 1600:// non-overlapping H2O band at 6.3 microns
    case 7250:// non-overlapping H2O band at 1.38 microns
	rhoalpha = rhoH2O*alpha;
	break;

    case 667: // non-overlapping CO2 band at 15  microns
    case 2410:// non-overlapping CO2 band at 4.3 microns
	rhoalpha = rhoCO2*alpha;
	break;

    case 3760:// H2O band at 2.7 microns, consisting of 3 overlapping sub-bands; with one overlapping CO2 band
    case 5350:// H2O+CO2 band; 2 overlapping bands: one H2O band at 1.87 microns and one CO2 band at 2.0 microns
	rhoalpha = rhoH2O*alpha + rhoCO2*overlapAlpha;
	break;

    default:
	WarningInFunction
	    << "  Unrecognized value of eta_c:" << eta_c
	    <<nl;
	break;
    }

    return rhoalpha;
}

/**************************************************************************\

\**************************************************************************/
void Foam::radiation::boxFixedBandsAbsorptionEmission::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aLambda
) const
{

    for (label j=0; j<nBands(); j++)
    {
	aLambda[j].primitiveFieldRef() = this->a(j);
    }
}

void Foam::radiation::boxFixedBandsAbsorptionEmission::correctEnFrac
(
    PtrList<volScalarField>& enFrac,
    blackBodyEmission& blackBody
) const
{
    enFrac[0] = 1;
    for (label j=1; j < nBands(); j++)
    {
	// radiating bands
        blackBody.correctEnFrac(enFrac[j], this->bands(j));
	// transparent composite band
	enFrac[0] -= enFrac[j];
    }
}
