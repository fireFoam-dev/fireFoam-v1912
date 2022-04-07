/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "absorptionCoeffs.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::radiation::absorptionCoeffs::~absorptionCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::absorptionCoeffs::checkT(const scalar T) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "using absorptionCoeffs out of temperature range:" << nl
            << "    " << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << nl << endl;
    }
}


const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::coeffs
(
    const scalar T
) const
{
    checkT(T);

    if (T < Tcommon_)
    {
        return lowACoeffs_;
    }
    else
    {
        return highACoeffs_;
    }
}


void Foam::radiation::absorptionCoeffs::initialise(const dictionary& dict)
{
    dict.readEntry("Tcommon", Tcommon_);
    dict.readEntry("Tlow", Tlow_);
    dict.readEntry("Thigh", Thigh_);
    dict.readEntry("invTemp", invTemp_);

    // luwi
    lowACoeffs_ = 0.;
    List<scalar> lowTdata(dict.lookup("loTcoeffs"));
    label ncoeffs = lowTdata.size();
    if( ncoeffs > nCoeffs_ )
    {
	FatalErrorIn
	(
	    "absorptionCoeffs::initialise(const dictionary& dict)"
	)   << "Number of coefficients ("<< ncoeffs <<") "
	    << "exceeds max. allowed ("<< nCoeffs_ <<")."
	    << exit(FatalError);
    }
    else
    {
	forAll(lowTdata,i)
	    lowACoeffs_[i] = lowTdata[i];
    }

    highACoeffs_ = 0.;
    List<scalar> highTdata(dict.lookup("hiTcoeffs"));
    ncoeffs = highTdata.size();
    if( ncoeffs > nCoeffs_ )
    {
	FatalErrorIn
	(
	    "absorptionCoeffs::initialise(const dictionary& dict)"
	)   << "Number of coefficients ("<< ncoeffs <<") "
	    << "exceeds max. allowed ("<< nCoeffs_ <<")."
	    << exit(FatalError);
    }
    else
    {
	forAll(highTdata,i)
	    highACoeffs_[i] = highTdata[i];
    }

}

// luwi
/* Parameters for boxFixedBandsAbsorptionEmission model, based on Exponential Wide Band model as described * \
\* in Section 11.10, table 11.3 of book: "Radiative Heat Transfer", 2013, by Michael F. Modest **************/
void Foam::radiation::absorptionCoeffs::initialiseEWB(const dictionary& dict)
{
    word species;
    dict.readEntry("species", species);
    dict.lookup("eta_c") >> eta_c_;
    label iauxError = 0;
    if (eta_c_==3760 && species=="H2O")
    {// allow specification of band strength parameter alpha_0 for up to 3 overlapping H2O bands at 3760
	List<scalar> alpha0list = dict.lookup("alpha_0");
	if ( alpha0list.size() > 3 )
	    iauxError = 3760;
	alpha_0_ = 0.;
	forAll(alpha0list, i)
	{
	    alpha_0_ += alpha0list[i];
	}
    }
    else
    {
	dict.lookup("alpha_0") >> alpha_0_;
    }

	
    // Sanity checks...   
    List<scalar> delta_k(dict.lookup("delta_k"));
    const dictionary overlap = dict.optionalSubDict("overlap");
    word overlapSpecies = "";
    scalar etac=-1.;

    if( overlap != dict )
    {
	overlap.lookup("species") >> overlapSpecies;
	overlap.lookup("eta_c") >> etac;
    }

    label ieta = int(eta_c_);
    
    switch(ieta)
    {
    case 7250:
	if ( species == "H2O" && delta_k == List<scalar>{1,0,1} && overlap == dict )
	    ieta = 0;
	break;
	
    case 5350:
	if ( species == "H2O" && delta_k == List<scalar>{0,1,1} )
	{//check for overlapping CO2 band
	    if ( overlap == dict )
	    {
		ieta = 0;
	    }
	    else if ( overlapSpecies=="CO2" && etac==5200 )
	    {
		overlap.lookup("alpha_0") >> overlapAlpha_0_;
		ieta = 0;
	    }
	}
	break;
	
    case 3760:
	if ( species == "H2O" && delta_k == List<scalar>{0,0,1} )
	{//check for overlapping CO2 band
	    if ( overlap == dict )
	    {
		ieta = 0;
	    }
	    else if ( overlapSpecies=="CO2" && etac==3660 )
	    {
		overlap.lookup("alpha_0") >> overlapAlpha_0_;
		ieta = 0;
	    }
	}
	break;

    case 2410:
	if ( species == "CO2" && delta_k == List<scalar>{0,0,1} && overlap == dict )
	    ieta = 0;
	break;

    case 1600:
	if ( species == "H2O" && delta_k == List<scalar>{0,1,0} && overlap == dict )
	    ieta = 0;
	break;

    case 667:
	if ( species == "CO2" && delta_k == List<scalar>{0,1,0} && overlap == dict )
	    ieta = 0;
	break;
	
    default:
	ieta = 9999;
	break;
    }

    // Print error messages
    const label& ierror = (iauxError? iauxError : ieta);
    switch (ierror)
    {
    case 0:
    {}
    break;
    
    case 7250:
	FatalErrorInFunction
	    << "  Band at "<<ieta<<" cm^-1 must be H2O band\n"
	    << "  with delta_k= (1 0 1) and no overlaps"
	    << exit(FatalError);	
	break;
	
    case 5350:
	FatalErrorInFunction
	    << "  Band at "<<ieta<<" cm^-1 must be H2O band\n"
	    << "  with delta_k= (0 1 1) and optional CO2 overlap at 5200 cm^-1"
	    << exit(FatalError);
	break;

    case 3760:
	FatalErrorInFunction
	    << "  Band at "<<ieta<<" cm^-1 must be H2O band\n"
	    << "  with delta_k= (0 0 1) and optional CO2 overlap at 3660 cm^-1\n"
	    << "  and up to 3 overlapping H2O sub-bands"
	    << exit(FatalError);
	break;

    case 2410:
	FatalErrorInFunction
	    << "  Band at "<<ieta<<" cm^-1 must be CO2 band\n"
	    << "  with delta_k= (0 0 1) and no overlaps"
	    << exit(FatalError);
	break;

    case 1600:
	FatalErrorInFunction
	    << "  Band at "<<ieta<<" cm^-1 must be H2O band\n"
	    << "  with delta_k= (0 1 0) and no overlaps"
	    << exit(FatalError);
	break;

    case 667:
	FatalErrorInFunction
	    << "  Band at "<<ieta<<" cm^-1 must be CO2 band\n"
	    << "  with delta_k= (1 0 1) and no overlaps"
	    << exit(FatalError);
	break;

    case 9999:
    default:
	FatalErrorInFunction
	    << "  Unrecognized band at "<<ieta<<" cm^-1\n"
	    << "  Available options are:\n"
	    << "  (\n"
	    << "    7250 (H2O)\n"
	    << "    5350 (H2O; with option of CO2 overlap at 5200 cm^-1)\n"
	    << "    3760 (H2O; with option of CO2 overlap at 3660 cm^-1)\n"
	    << "    2410 (CO2)\n"
	    << "    1600 (H2O)\n"
	    << "     667 (CO2)\n"
	    << "  )\n"
	    << exit(FatalError);
	break;
    }
}

// ************************************************************************* //
