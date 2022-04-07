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

#include "wsggAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggAbsorptionEmission,
            dictionary
        );
    }
}

namespace Foam
{
    const Foam::Enum
    <
        Foam::radiation::wsggAbsorptionEmission::correlationType
    >correlationTypeNames_
    {
        { Foam::radiation::wsggAbsorptionEmission::correlationType::absCoeffJohansson, "johansson" },
        { Foam::radiation::wsggAbsorptionEmission::correlationType::absCoeffFranca, "franca" },
        { Foam::radiation::wsggAbsorptionEmission::correlationType::absCoeffSmith, "smith" }
    };
}

//
const Foam::Switch
Foam::radiation::wsggAbsorptionEmission::getGrayFlag(const dictionary& dict) const
{
    const word userModel = dict.lookupOrDefault<word>("spectralMethod","gray");
    Switch grayFlag = true;

    if( userModel == "gray" || userModel == "grey" )
	grayFlag = true;
    else if( userModel == "banded" )
	grayFlag = false;
    else
    {
	FatalErrorIn
	(
	    "radiation::wsggAbsorptionEmission::checkGrayModel(const dictionary& dict)"
        )   << "Unrecognized spectralMethod: "<< userModel << nl
	    << "Valid options are:"<<nl
	    <<"("<<nl
	    <<"   'gray' (or 'grey') "<<nl
	    <<"   'banded' " << nl
	    <<")"
	    << exit(FatalError);
    }

    return grayFlag;
}

const Foam::dictionary
Foam::radiation::wsggAbsorptionEmission::wsggDict(const fvMesh& mesh)
{
    // Read wsgg correlation parameters from lookup file
    IOdictionary wsggFileDict
    (
	IOobject
	(
	    wsggFileName_,
	    mesh.time().constant(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	)
    );

    if( correlationID_=="" )
    {
	forAllConstIter(dictionary, wsggFileDict.subDict( wsggModelName_ ), iter)
	{
	    if(iter().isDict())
	    {
		correlationID_ = iter().keyword();
		break;
	    }
	}
    }

    return( (wsggFileDict.subDict( wsggModelName_ )).subDict( correlationID_ ) );

}

const Foam::psiReactionThermo&
Foam::radiation::wsggAbsorptionEmission::thermoRef(const fvMesh& mesh) const
{
    return ( mesh.lookupObject<psiReactionThermo>("thermophysicalProperties") );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggAbsorptionEmission::wsggAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    wsggModelName_(coeffsDict_.lookupOrDefault<word>("correlationType","johansson")),
    correlationID_(coeffsDict_.lookupOrDefault<word>("correlationID","") ),
    wsggModel_(	correlationTypeNames_[ wsggModelName_ ] ),
    wsggFileName_( coeffsDict_.lookupOrDefault<word>("lookupFile","wsggProperties") ),
    wsggDict_( wsggDict(mesh) ),
    gasThermo_(	dynamic_cast<const reactingMixture<gasHThermoPhysics>&>	( thermoRef(mesh) ) ),
    p_( thermoRef(mesh).p() ),
    T_( thermoRef(mesh).T() ),
    Y_( thermoRef(mesh).composition().Y() ),
    indexCO2_( gasThermo_.species()["CO2"] ),
    indexH2O_( gasThermo_.species()["H2O"] ),
    molWtCO2_( gasThermo_.speciesData()[indexCO2_].W() ),
    molWtH2O_( gasThermo_.speciesData()[indexH2O_].W() ),
    Csoot_(coeffsDict_.lookupOrDefault<scalar>("Csoot",0.0)),
    mixtureFractionSoot_( Csoot_==0? false : true ),
    isGray_( getGrayFlag(coeffsDict_) ),
    pathLength_
    (
	coeffsDict_.lookupOrDefault<scalar>("pathLength",1.0)
	*dimensionedScalar("unity", dimLength, 1.0)
    ),
    Tref_
    (
	wsggDict_.lookupOrDefault<scalar>("Tref",1.0)
       *dimensionedScalar("one",dimTemperature,1.0)
    )

{

    label nBand = 1; //Always includes one clear gas
    forAllConstIter(dictionary, wsggDict_, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }

	const word& key = iter().keyword();

	for(label iband=0; iband<maxBands_; iband++)
	{
	    if(
		  ( key == "greyGas"+name(iband) )
	       || ( key == "grayGas"+name(iband) )
	      )
	    {
		const dictionary& bandDict = iter().dict();
		coeffs_[nBand].initialise(bandDict);
		nBand++;
		break;
	    }
	}
    }
    nBands_ = nBand;

    if( nBands_ == 1 )
    {
	WarningInFunction << "No gray gas specified! \n"
			  << "\t Will assume transparent gas mixture."
			  << nl;
    }
    else
    {
	Info<<"WSGG correlation:"<<nl
	    <<"\t" <<wsggModelName_ <<"--" << correlationID_
	    <<" with " <<(nBands_-1)
	    <<" gray gases and 1 transparent gas."
	    <<nl;
    }

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggAbsorptionEmission::~wsggAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggAbsorptionEmission::aCont(const label bandI) const
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
    static volScalarField mixMolWt(ta.ref()*dimensionedScalar("one", dimLength, 1.0));
    static volScalarField X_CO2(mixMolWt);
    static volScalarField X_H2O(mixMolWt);

    if( bandI==0 )
    {// Initialization
	mixMolWt = 0.*ggWeight(bandI);
	forAll(Y_,specieI)
	{
	    const scalar& molWti = gasThermo_.speciesData()[specieI].W();
	    mixMolWt += Y_[specieI]/molWti;
	}
	mixMolWt = 1/mixMolWt;

	X_CO2 = mixMolWt*Y_[indexCO2_]/molWtCO2_;
	X_H2O = mixMolWt*Y_[indexH2O_]/molWtH2O_;

	//Gray calculation
	if( isGray_ )
	{
	    volScalarField emissivity(0.0*X_CO2);
	    for(label iband=1; iband<nBands_; iband++)
	    {
		volScalarField ai("ai", aCont(iband));
		emissivity += ggWeight(iband) * ( 1. - exp(-pathLength_ * ai) );
	    }
	    emissivity.min(0.9999);
	    a = -log( 1. - emissivity )/pathLength_;
	    //-- add soot
	    if( mixtureFractionSoot_ )
	    {
		volScalarField ggEnFrac("ggEnFrac", (X_CO2-X_CO2+1.0));
		for(label iband=1; iband<nBands_; iband++)
		{
		    const volScalarField& enFraci = ggWeight(iband);
		    ggEnFrac -= enFraci;
		}

		const volScalarField& fv =mesh_.lookupObject<volScalarField>("soot");
		a += ggEnFrac * Csoot_ * fv * T_
		    /dimensionedScalar("unity",dimTemperature/dimLength,1.0);
	    }
	}
    }
    else
    {// Gray gas absorption coefficients
	switch( wsggModel_ )
	{
	   case absCoeffJohansson:
	   default:
	   {
	       //b[0] in units of 1/(bar*m); where bar[=] sum of partial pressures of H2O+CO2
	       const scalar K1 = coeffs_[bandI].coeffs(T_[0])[0];
	       const scalar K2 = coeffs_[bandI].coeffs(T_[0])[1];
	       volScalarField MR("MR", max( X_H2O, VSMALL ) / max( X_CO2, SMALL ));
	       a = (K1 + K2*MR) * (X_H2O + X_CO2) * p_ * paToBar_;
	   }
	   break;

	   case absCoeffFranca:
	   case absCoeffSmith:
	   {
	       //k in units of 1/(atm*m); where atm[=] sum of partial pressures of H2O+CO2
	       //kappa is a constant for each gray gas
	       const scalar k = (coeffs_[bandI].coeffs(T_[0]))[0];
	       a = k*(X_H2O+X_CO2)*p_*paToAtm_;
	   }
	   break;
	}

	if( mixtureFractionSoot_ && !isGray_ )
	{
	    const volScalarField& fv =mesh_.lookupObject<volScalarField>("soot");
	    a += ggWeight(bandI) * Csoot_ * fv * T_
		/dimensionedScalar("unity",dimTemperature/dimLength,1.0);
	}
    }


    return ta;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggAbsorptionEmission::ggWeight(const label bandI) const
{

    tmp<volScalarField> tggWeight
    (
        new volScalarField
        (
            IOobject
            (
                "ggWeight",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("ggWeight", dimless, 1.0)
        )
    );

    static volScalarField mixMolWt(tggWeight.ref());
    static volScalarField X_CO2(mixMolWt);
    static volScalarField X_H2O(mixMolWt);
    static volScalarField Tdimless(mixMolWt);

    if( bandI==0 )
    {// Initialization
	mixMolWt = 0.;
	forAll(Y_,specieI)
	{
	    const scalar& molWti = gasThermo_.speciesData()[specieI].W();
	    mixMolWt += Y_[specieI]/molWti;
	}
	mixMolWt = 1/mixMolWt;

	X_CO2 = mixMolWt*Y_[indexCO2_]/molWtCO2_;
	X_H2O = mixMolWt*Y_[indexH2O_]/molWtH2O_;
	Tdimless = T_/Tref_;
    }
    else
    {// Gray gas weights
	scalarField& enFrac = tggWeight.ref().primitiveFieldRef();
	//--Johansson, Franca, Smith coeffs b[] are independent of T
	const absorptionCoeffs::coeffArray& b = coeffs_[bandI].coeffs(T_[0]);

	switch( wsggModel_ )
	{
	   case absCoeffJohansson:
	   default:
	   {
	       volScalarField MR("MR", max( X_H2O, VSMALL ) / max( X_CO2, SMALL ));
	       volScalarField MR2("MR2", MR*MR);
	       const volScalarField c1("c1", b[2] + b[5]*MR + b[ 8]*MR2);
	       const volScalarField c2("c2", b[3] + b[6]*MR + b[ 9]*MR2);
	       const volScalarField c3("c3", b[4] + b[7]*MR + b[10]*MR2);

	       enFrac = (c3*Tdimless + c2)*Tdimless + c1;
	   }
	   break;

	   case absCoeffSmith:
	   {
	       enFrac = ((b[4]*Tdimless + b[3])*Tdimless + b[2])*Tdimless + b[1];
	   }
	   break;

	   case absCoeffFranca:
	   {
	       enFrac = (((b[5]*Tdimless + b[4])*Tdimless + b[3])*Tdimless + b[2])*Tdimless + b[1];
	   }
	   break;
	}

    }

    return tggWeight;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggAbsorptionEmission::eCont(const label bandI) const
{
    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggAbsorptionEmission::ECont(const label bandI) const
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

void Foam::radiation::wsggAbsorptionEmission::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aLambda

) const
{
    a = dimensionedScalar("zero", dimless/dimLength, 0.0);

    for(label j=0; j<nBands(); j++)
	aLambda[j].primitiveFieldRef() = this->a(j);

}

void Foam::radiation::wsggAbsorptionEmission::correctEnFrac
(
    PtrList<volScalarField>& enFrac,
    blackBodyEmission& blackBody
) const
{
    enFrac[0] = ggWeight(0);

    if( !isGrey() )
    {
	for (label j=1; j < nBands_; j++)
	{
	    enFrac[j] = this->ggWeight(j);
	    enFrac[0] -= enFrac[j];
	}
    }

}




// ************************************************************************* //
