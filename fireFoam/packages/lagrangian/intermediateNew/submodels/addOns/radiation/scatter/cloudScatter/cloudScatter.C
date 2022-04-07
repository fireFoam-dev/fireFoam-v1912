/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "cloudScatter.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(cloudScatter, 0);

        addToRunTimeSelectionTable
        (
            scatterModel,
            cloudScatter,
            dictionary
        );
    }
}

namespace Foam
{
    Foam::Enum
    <
        Foam::radiation::cloudScatter::phaseFunctionType
    >
    Foam::radiation::cloudScatter::phaseFunctionTypeNames_
    {
        {phaseFunctionType::preCorrected, "preCorrected"},
        {phaseFunctionType::henyeyGreenstein, "henyeyGreenstein"},
        {phaseFunctionType::hgDirect, "hgdirect"},
	{phaseFunctionType::hgfixedAsymmetry, "hgfixedAsymmetry"}
    };
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::cloudScatter::cloudScatter
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    scatterModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    cloudNames_(coeffsDict_.lookup("cloudNames")),
//<fmglobal>
    phaseFoption_
    (
        phaseFunctionTypeNames_[coeffsDict_.lookupOrDefault<word>("phaseFunction","preCorrected")]
    ),
    asymmetryFactor_
    (
	coeffsDict_.lookupOrDefault<scalar>("asymmetryFactor",0.9)
    ),
    nminTheta_( coeffsDict_.lookupOrDefault<label>("nminTheta",32) ),
    nminPhi_( coeffsDict_.lookupOrDefault<label>("nminPhi",16) ),
    nphaseFUniq_(0),
    ncscatUniq_(0),
    doInit_(true)
//</fmglobal>
{
//<fmglobal>
    if ( cloudNames_.size() > 1 )
    {
    FatalErrorIn
        (
        "radiation::cloudScatter()"
        )   << "cloudScatter model will only work for 1 cloud "
                << "(" << cloudNames_.size() << " clouds specified)."
                << exit(FatalError);
    }
//</fmglobal>
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::cloudScatter::~cloudScatter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::radiation::cloudScatter::sigmaEff(const label bandI, const label iRay, const scalar omega) const //ankur
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
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

    forAll(cloudNames_, i)
    {
        const thermoCloud& tc
        (
            mesh_.objectRegistry::lookupObject<thermoCloud>(cloudNames_[i])
        );
        // luwi
        if( tc.radiation() )
        {
	    if (iRay<0)
	    {
		tsigma.ref() = tc.sigmap(bandI);
		return tsigma;
	    }
	    else
	    {
		switch( phaseFoption_ )
		{
		case preCorrected:
		default:
		{//Tabulated scattering coefficients pre-corrected for self-in-scattering
		    tsigma.ref() = tc.sigmap(bandI); // ankur
		}
		break;

		case hgDirect:
		{
		    tsigma.ref() = tc.sigmap(bandI).ref() *
			(
			    1.
			    - omega*tc.phaseFuncHGDirect(bandI,iRay,iRay)
			    /(4.*constant::mathematical::pi)
			    );
		}
		break;

		case henyeyGreenstein:
		{//No pre-correction; correct for self-in-scattering here
		    const label& icscatRow = cscatRow_[iRay];
		    const scalarField& cscat = cscatUniq_[bandI][icscatRow].field();

		    tsigma.ref().primitiveFieldRef() = tc.sigmap(bandI).ref().primitiveFieldRef()*cscat;

		}
		break;

		case hgfixedAsymmetry:
		{//No pre-correction; correct for self-in-scattering here
		    const label& icscatRow = cscatRow_[iRay];
		    const scalar& cscat = refCscatUniq_[icscatRow];
		    tsigma.ref().primitiveFieldRef() = tc.sigmap(bandI).ref().primitiveFieldRef()*cscat;
		}
		break;

		}
	    }
	}
    }

    // return 3.0*tsigma;
    // ankur, removing 3.0 here, to make it generic for P1 and fvDOM
    return tsigma; // luwi
}

//<fmglobal>
Foam::tmp<Foam::volScalarField>
Foam::radiation::cloudScatter::pFunc(const label bandI, const label sour, const label dest) const
{
    tmp<volScalarField> tPhaseFunc
    (
        new volScalarField
        (
            IOobject
            (
                "phaseFunc",
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

    forAll(cloudNames_, i)
    {
        const thermoCloud& tc
            (
             mesh_.objectRegistry::lookupObject<thermoCloud>(cloudNames_[i])
            );
        // Compute only if radiation was activated for cloud
        if( tc.radiation() )
        {
            switch(phaseFoption_)
            {
                case preCorrected:
                default:
                    tPhaseFunc.ref() = tc.sigmap(bandI)*tc.phaseFunc(bandI,sour,dest);   // sigma_cloud * phaseFunc_cloud [sigma not to be included
                    // while solving for RTE since it is already included here]
                    break;
                case hgDirect:
                    {
                        if( sour != dest )
                            tPhaseFunc.ref() = tc.sigmap(bandI)*tc.phaseFuncHGDirect(bandI,sour,dest);
                    }
                    break;
                case henyeyGreenstein:
                    {
                        if( sour != dest )
                        {
                            const label& iuniqF = phaseFindx_[sour][dest];

			    tPhaseFunc.ref().primitiveFieldRef() = phaseFUniq_[bandI][iuniqF].field();
                        }
                    }
                    break;
	        case hgfixedAsymmetry:
		    {
			if ( sour != dest )
			{
			    const label& iuniqF = phaseFindx_[sour][dest];
			    tPhaseFunc.ref().primitiveFieldRef() = refPhaseFUniq_[iuniqF];
			}
		    }
		    break;
            }
        }
    }

    return tPhaseFunc;

}

void
Foam::radiation::cloudScatter::initPhaseFuncs(const fvDOM& dom)
{
    nRays_ = dom.nRay();
    label angleAverage =
    (
	phaseFoption_ == hgfixedAsymmetry
     &&	mesh_.nSolutionD() == 3
     && dom.nTheta() < nminTheta_
     && dom.nPhi() < nminPhi_
    );


    if (angleAverage)
    {
	angleAveragedPhaseF(dom);
    }
    else
    {
	setupMaps(dom);
    }

    if (phaseFoption_ != hgfixedAsymmetry)
    {
	setupMaps(dom);
	nBands_ = dom.nLambda();
	cscatUniq_.setSize(nBands_);
	phaseFUniq_.setSize(nBands_);
	const label& nphaseFUniq(nphaseFUniq_);
	const label& ncscatUniq(ncscatUniq_);

	forAll(cscatUniq_, lambdaI)
	{
	    cscatUniq_.set(lambdaI, new PtrList<volScalarField>(ncscatUniq));
	    phaseFUniq_.set(lambdaI, new PtrList<volScalarField>(nphaseFUniq));
	    forAll(cscatUniq_[lambdaI], jj)
	    {
		cscatUniq_[lambdaI].set
		(
		    jj,
		    new volScalarField
		    (
			IOobject
			(
			    "cscatUniq_" + Foam::name(lambdaI) + "_" + Foam::name(jj),
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
	    forAll(phaseFUniq_[lambdaI], jj)
	    {
		phaseFUniq_[lambdaI].set
		(
		    jj,
		    new volScalarField
		    (
			IOobject
			(
			    "phaseFUniq_" + Foam::name(lambdaI) + "_" + Foam::name(jj),
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
	}
    }

}

void Foam::radiation::cloudScatter::angleAveragedPhaseF(const fvDOM& dom)
{
    // Computes average phase function for radiative scattering between two solid angles,
    // using a finer discretization than the one solved in fvDOM.
    // The finer resolution is used only for averaging; the results are mapped
    // back to the fvDOM discretization.
    // This is done only once; phase function values remain fixed throughout simulation.

    //first, determine resolution of fine discretization
    label ndivTheta = ceil( scalar(nminTheta_)/dom.nTheta() );
    label ndivPhi = ceil( scalar(nminPhi_)/dom.nPhi() );

    label nTheta = dom.nTheta()*ndivTheta;
    label nPhi = dom.nPhi()*ndivPhi;

    label nRays = 4*nPhi*nTheta;
    const scalar& pi(constant::mathematical::pi);
    scalar deltaPhi = pi/(2.0*nPhi);
    scalar deltaTheta = pi/nTheta;

    // Setup fine angular discretization for integration (need dAve and omega)
    List<vector> dAve(nRays);
    scalarList omega(nRays);
    label i = 0;
    for (label n = 1; n <= nTheta; n++)
    {
	scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
	scalar sinTheta = sin(thetai);

	for (label m=1; m <= 4*nPhi; m++)
	{
	    scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
	    omega[i] = 2.0*sinTheta*sin(deltaTheta/2.0)*deltaPhi;

	    scalar sinPhi = sin(phii);
	    scalar cosPhi = cos(phii);

	    dAve[i] = vector
	    (
		sinPhi
		*Foam::sin(0.5*deltaPhi)
		*(deltaTheta - Foam::cos(2.0*thetai)
		  *Foam::sin(deltaTheta)),
		cosPhi
		*Foam::sin(0.5*deltaPhi)
		*(deltaTheta - Foam::cos(2.0*thetai)
		  *Foam::sin(deltaTheta)),
		0.5*deltaPhi*Foam::sin(2.0*thetai)*Foam::sin(deltaTheta)
	    );
	    dAve[i] /= mag(dAve[i]);

	    i++;
	}
    }

    // Compute phase function for fine discretization...
    label nmax = (nRays*(nRays-1))/2 + 1;//symmetric matrix, with all diagonal elements identical
    scalarList phaseFUniq(nmax,0.);
    List<labelList> phaseFindx(nRays,labelList(nRays,0));
    const scalar& g(asymmetryFactor_);// sample asymmetry factor value
    const scalar myzero(1.e-6);

    label nphaseFUniq=0;
    for(label iRow=0; iRow<nRays; iRow++)
    {
	for(label iCol=0; iCol<nRays; iCol++)
	{
	    if(iCol<iRow)
	    {//matrix lower triangle
		phaseFindx[iRow][iCol] = phaseFindx[iCol][iRow];
	    }
	    else if ( iCol>iRow )
	    {//matrix upper triangle, excluding diagonal
		scalar cosTheta = dAve[iRow] & dAve[iCol];
		scalar rval0 = (1. - g*g)/(Foam::pow((1. + g*g - 2.*g*cosTheta),1.5));
		label iuniq = 1;
		for(label jj=0; jj<nphaseFUniq; jj++)
		{
		    if( fabs(rval0 - phaseFUniq[jj])< myzero )
		    {
			phaseFindx[iRow][iCol] = jj;
			iuniq=0;
			break;
		    }
		}
		if( iuniq )
		{
		    phaseFUniq[nphaseFUniq] = rval0;
		    phaseFindx[iRow][iCol] = nphaseFUniq;
		    nphaseFUniq++;
		}
	    }
	}
    }

    // Then, do integral to compute solid-angle-averaged phaseF
    // ...first, map local (for phaseF solid-angle-averaging) discretization to fvDOM discretization
    label iRay=0;
    label dom_iRay=0;
    label dom_n=0;
    labelList rayMap(nRays, 0);
    label ival0=0;
    for (label n = 0; n < nTheta; n++)
    {
	if( ! (n % ndivTheta ) )
	{
	    ival0 = dom_n*(4*dom.nPhi());
	    dom_n++;
	}

	label dom_m = 0;
	for (label m=0; m < 4*nPhi; m++)
	{
	    if( ! (m % ndivPhi) )
	    {
		dom_iRay = ival0 + dom_m;
		dom_m++;
	    }
	    rayMap[iRay] = dom_iRay;

	    iRay++;
	}
    }

    // ...setup local storage for phaseF in fvDOM discretization
    List<scalarList> dom_phaseF(nRays_,scalarList(nRays_,0.));
    // ...initialize class arrays
    nmax = (nRays_*(nRays_-1))/2 + 1;
    refPhaseFUniq_.setSize(nmax);
    label ncscat = max(1,ceil(dom.nTheta()/2.));
    refCscatUniq_.setSize(ncscat);
    cscatRow_.setSize(nRays_);
    phaseFindx_.setSize(nRays_, cscatRow_);
    // ...initialize counters
    nphaseFUniq_ = ncscatUniq_ = 0;

    // ...compute double integral {phaseF dOmega' dOmega}
    for (iRay=0; iRay<nRays; iRay++)
    {
	label iRow = rayMap[iRay];
	for (label jRay=0; jRay<nRays; jRay++)
	{
	    label iCol = rayMap[jRay];
	    if ( iCol > iRow )
	    {// matrix upper triangle
		dom_phaseF[iRow][iCol] += phaseFUniq[ phaseFindx[iRay][jRay] ] * omega[iRay] * omega[jRay];
	    }
	}
    }

    // ...and finally, compute solid-angle-averages in fvDOM discretization
    for (label iRow=0; iRow<nRays_; iRow++)
    {
	scalar rval = 0.0;
	scalar iomega = dom.IRay(iRow).omega();
	for (label iCol=0; iCol<nRays_; iCol++)
	{
	    if(iCol<iRow)
	    {//matrix lower triangle
		phaseFindx_[iRow][iCol] = phaseFindx_[iCol][iRow];
		dom_phaseF[iRow][iCol] = dom_phaseF[iCol][iRow];
	    }
	    else if ( iCol>iRow )
	    {//matrix upper triangle
		dom_phaseF[iRow][iCol] /= (iomega*dom.IRay(iCol).omega());
		const scalar& rval0( dom_phaseF[iRow][iCol] );
		label iuniq=1;
		for(label jj=0; jj<nphaseFUniq_; jj++)
		{
		    if( fabs(rval0 - refPhaseFUniq_[jj])< myzero )
		    {
			phaseFindx_[iRow][iCol] = jj;
			iuniq=0;
			break;
		    }
		}
		if( iuniq )
		{
		    refPhaseFUniq_[nphaseFUniq_] = rval0;
		    phaseFindx_[iRow][iCol] = nphaseFUniq_;
		    nphaseFUniq_++;
		}
	    }
	    //compute diagonal to enforce conservation of scattered energy
	    rval += dom_phaseF[iRow][iCol]*dom.IRay(iCol).omega();
	}
	rval /= (4*pi);

	label iuniq=1;
	for(label jj=0; jj<ncscatUniq_; jj++)
	{
	    if (fabs(rval - refCscatUniq_[jj]) < myzero )
	    {
		cscatRow_[iRow] = jj;
		iuniq = 0;
		break;
	    }
	}
	if(iuniq)
	{
	    refCscatUniq_[ncscatUniq_] = rval;
	    cscatRow_[iRow] = ncscatUniq_;
	    ncscatUniq_++;
	}
    }

    //resize arrays
    refPhaseFUniq_.setSize(nphaseFUniq_);
    refCscatUniq_.setSize(ncscatUniq_);
}

void Foam::radiation::cloudScatter::setupMaps(const fvDOM& dom)
{
//- Setup maps for efficient computation of in-scattering
//  Exploit matrix symmetry (Phi_i,j = Phi_j,i), plus other features of the
//  angular discretization that yield many repeated values of the matrix
//  Find where the unique values are, as well as the matrix row and column
//  indices i and j that share these values
    label nmax = (nRays_*(nRays_-1))/2 + 1;//symmetric matrix, with all diagonal elements identical
    scalarList& phaseFUniq( refPhaseFUniq_ );
    phaseFUniq.setSize(nmax);
    scalarList cosTheta(nmax,0.);
    scalarList magdAve(nRays_, 0.);
    forAll(magdAve, i)
    {
	magdAve[i] = mag(dom.IRay(i).dAve());
    }
    List<labelList> iphaseFUniq(nmax,labelList(label(2),0));
    label nphaseFUniq=0;
    label ncscatUniq(0);
    label ncscat = max(1,ceil(dom.nTheta()/2.));
    labelList icscatUniq(ncscat, 0);
    scalarList& cscatUniq( refCscatUniq_ );
    cscatUniq.setSize(ncscat);
    scalar g = asymmetryFactor_;// sample asymmetry factor value
    const scalar myzero(1.e-6);
    cscatRow_.setSize(nRays_);
    phaseFindx_.setSize(nRays_, cscatRow_);
    for(label iRow=0; iRow<nRays_; iRow++)
    {
    scalar rval(0.);
    for(label iCol=0; iCol<nRays_; iCol++)
    {
	// skip diagonal
	if(iCol==iRow)
	    continue;

        if(iCol<iRow)
        {//matrix lower triangle
        phaseFindx_[iRow][iCol] = phaseFindx_[iCol][iRow];
        }
        else
        {//matrix upper triangle, including diagonal
	scalar rtemp = dom.IRay(iRow).dAve() & dom.IRay(iCol).dAve();
	rtemp /= ( magdAve[iRow] * magdAve[iCol] );
        scalar rval0 = (1. - g*g)/(Foam::pow((1. + g*g - 2.*g*rtemp),1.5));
        label iuniq=1;
        for(label jj=0; jj<nphaseFUniq; jj++)
        {
            if( fabs(rval0 - phaseFUniq[jj])< myzero )
            {
            phaseFindx_[iRow][iCol] = jj;
            iuniq=0;
            break;
            }
        }
        if( iuniq )
        {
	    cosTheta[nphaseFUniq] = rtemp;
            phaseFUniq[nphaseFUniq] = rval0;
            phaseFindx_[iRow][iCol] = nphaseFUniq;
            iphaseFUniq[nphaseFUniq][0] = iRow;
            iphaseFUniq[nphaseFUniq][1] = iCol;
            nphaseFUniq++;
        }
        }

        label iind=phaseFindx_[iRow][iCol];
        rval += phaseFUniq[iind]*dom.IRay(iCol).omega();
    }

    rval *= ( 1.0/(4.*constant::mathematical::pi) );
    label iuniq=1;
    for(label jj=0; jj<ncscatUniq; jj++)
    {
        if (fabs(rval - cscatUniq[jj]) < myzero )
        {
        cscatRow_[iRow] = jj;
        iuniq = 0;
        break;
        }
    }
    if(iuniq)
    {
        cscatUniq[ncscatUniq] = rval;
        cscatRow_[iRow] = ncscatUniq;
        icscatUniq[ncscatUniq] = iRow;
        ncscatUniq++;
    }
    }

    nphaseFUniq_ = nphaseFUniq;
    // resize refPhaseFUniq_
    refPhaseFUniq_.setSize(nphaseFUniq_);

    ncscatUniq_ = ncscatUniq;
    icscatUniq_.setSize(ncscatUniq_);
    refCscatUniq_.setSize(ncscatUniq_); //resized
    forAll(icscatUniq_, i)
    icscatUniq_[i] = icscatUniq[i];
    cosTHETA_.setSize(nphaseFUniq);
    forAll(cosTHETA_, i)
	cosTHETA_[i] = cosTheta[i];

    iphaseFUniqRow_.setSize(nphaseFUniq);
    iphaseFUniqCol_.setSize(nphaseFUniq);
    forAll(iphaseFUniqRow_, i)
    {
    iphaseFUniqRow_[i] = iphaseFUniq[i][0];
    iphaseFUniqCol_[i] = iphaseFUniq[i][1];
    }

}

void
Foam::radiation::cloudScatter::updatePhaseFuncs(const fvDOM& dom)
{
    //First select the thermo cloud;
    //NOTE: this implementation works only for simulations with a single cloud
    const thermoCloud& tc( mesh_.objectRegistry::lookupObject<thermoCloud>(cloudNames_[0]) );
    //Return if radiation is not active for cloud; or if not using henyeyGreenstein
    if ( !tc.radiation() || (phaseFoption_ == preCorrected) || (phaseFoption_ == hgDirect) )
    return;

    if ( doInit_ )
    {
	initPhaseFuncs(dom);
	doInit_ = 0;
    }
    // Nothing to compute for hgfixedAsymmetry...
    if ( phaseFoption_==hgfixedAsymmetry )
    {
	return;
    }

    tmp<volScalarField> asyFac
    (
        new volScalarField
        (
            IOobject
            (
                "asyFac",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& asyFacVal = asyFac.ref().primitiveFieldRef();

    //   Compute and store all phase functions phaseF_ij

    for(label bandI=0; bandI<nBands_; bandI++)
    {
    asyFac.ref() = tc.asyFac(bandI);
    scalarField asyFacVal2(asyFacVal*asyFacVal);
    //- First, compute all unique phaseF_ij
    for (label iuniqF=0; iuniqF < nphaseFUniq_; iuniqF++)
    {
        const scalar& cosTHETA = cosTHETA_[iuniqF];

        phaseFUniq_[bandI][iuniqF].field() =
        (1. - asyFacVal2)/Foam::pow((1. + asyFacVal2 - 2.*asyFacVal*cosTHETA),1.5);

    }
    //- Then, compute normalization factors
    for(label iuniq=0; iuniq<ncscatUniq_; iuniq++)
    {
        const label& iRow = icscatUniq_[iuniq];
        scalarField& cscatUniq = cscatUniq_[bandI][iuniq].field();
        cscatUniq = 0.0;

        for(label iCol=0; iCol<nRays_; iCol++)
        {
        //skip diagonal
	if(iCol==iRow)
	    continue;

        const label& iuniqF = phaseFindx_[iRow][iCol];

        cscatUniq += phaseFUniq_[bandI][iuniqF].field()*dom.IRay(iCol).omega();
        }
        cscatUniq *=
        (
            1.0
            / (4.*constant::mathematical::pi)
        );
    }
    }

}

//Flag to indicate cloud model is being used with preCorrected
Foam::Switch
Foam::radiation::cloudScatter::cloudCorrectedSigma() const
{
    return ( phaseFoption_ == preCorrected );
}
//</fmglobal>

// ************************************************************************* //
