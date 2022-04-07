// IMPORTANT: Used in computations (not just diagnostics)...
if ( !skipInfoOutputSolverCalcs )
{
    //- Convective heat flux
    convectiveHeatFlux_T = 
	- fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);        

	convectiveHeatFlux_L = 
	    - fvc::interpolate(turbulence->alpha())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);
	
    //- Solid mass loss rate
    forAll(QcWallFunction.boundaryField(), patchi)
    {
	    if (mesh.boundary()[patchi].type() == "mappedWall") 
	    {
                mlrSolid.boundaryFieldRef()[patchi] =
                (
                    - phi.boundaryField()[patchi]
                     /mesh.boundary()[patchi].magSf()
                     *1000.
                );
			QcWallFunction.boundaryFieldRef()[patchi] = convectiveHeatFlux_T.boundaryField()[patchi];
	    }
    }
    
     //- Flame Extinction (Bert)
     volFracSpray = parcels.theta();
     rhoSpray = rhoWater*volFracSpray;
}
   //- For sprinkler skipping
   sprayFluxNet = parcels.sprayFluxV();
   sprayFluxAF  = parcels.sprayFluxS();
   sprayTemperature = parcels.sprayTemperature();

   //- Needed for lattice melt
   sprayMass = parcels.sprayMass();

   //- Update spray BC
   sprayFluxNet.correctBoundaryConditions(); 
   sprayMass.correctBoundaryConditions();
   sprayTemperature.correctBoundaryConditions();

   //- For flame visualization
   if (thermo.composition().contains("O2")) 
   {
       ft = (fu*s-O2+YO2Inf)/(s+YO2Inf);
   }
   else
   {
       ft = fu;
   }

   ft.max(0.0);
   ft.min(1.0);
    

/////////////////////////////////////////////////////////////////////////////////////////////////

// For diagnostics only...
if
(
    writeInfoOutput 
 && !(runTime.timeIndex() % infoOutputInterval)
)
{
    // Sometimes used, for diagnostics...
    tmp<fv::convectionScheme<scalar> > mvConvection
    (
	fv::convectionScheme<scalar>::New
	(
	    mesh,
	    fields,
	    phi,
	    mesh.divScheme("div(phi,Yi_h)")
        )
    );
        
    HRR_fu = ( - mvConvection->interpolate(phi,fu)*phi
	       + (
		   fvc::interpolate
		   (
		       turbulence->alphaEff()*fvc::grad(fu)
                   )
		   & mesh.Sf()
                 )                         
             ) * qFuel;

    
    phiFt = mvConvection->interpolate(phi,ft)*phi
	    - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(ft)*mesh.magSf();
    
    phiO2 = mvConvection->interpolate(phi,O2)*phi
            - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(O2)*mesh.magSf();

    phiCO2 = mvConvection->interpolate(phi,CO2)*phi
            - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(CO2)*mesh.magSf();

    hc = thermo.hc();
    phiHc = phi*mvConvection->interpolate(phi,hc)
            - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(hc)*mesh.magSf();

    phiHs = phi*mvConvection->interpolate(phi,thermo.he())
            - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(thermo.he())*mesh.magSf();

    h = thermo.he() + hc;
    phiH  = phi*mvConvection->interpolate(phi,h)
            - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(h)*mesh.magSf();

    // Rarely used, for diagnostics...
    if ( verboseInfoOutput )
    {
	momentumX = phi * fvc::interpolate(U)->component(0); // * runTime.deltaTValue() * unitMom; 
	momentumY = phi * fvc::interpolate(U)->component(1); // * runTime.deltaTValue() * unitMom; 
	momentumZ = phi * fvc::interpolate(U)->component(2); // * runTime.deltaTValue() * unitMom; 	
	
	UT = U * T;
	rhoU = rho * U;
	rhoT = rho * T;
	rhoFu = rho * fu; 
	
	B = turbulence->R();
	
	if (!constD) 
	{
	    d = turbulence->alpha()/lewisNo/rho;
	}
	else 
	{
	    d = DM;
	}
	dSgs = (turbulence->alphaEff()-turbulence->alpha())/rho;

	//forAll(QcWallFunction.boundaryField(), patchi)
	//{
	//    if (mesh.boundary()[patchi].type() == "mappedWall") 
	//    {
	//	forAll ( QcWallFunction.boundaryField()[patchi],faceI)
	//	{
	//	    scalar mlr = mlrSolid.boundaryFieldRef()[patchi][faceI] * 2.5;

	//	    if (mlr < 0.1)
	//	    {
	//		QcWallFunction.boundaryFieldRef()[patchi][faceI] = min(max(0,convectiveHeatFlux_L.boundaryField()[patchi][faceI]),QcThreshold)/QcThreshold*QcFlame;
	//	    }
	//	    else
	//	    {
	//		QcWallFunction.boundaryFieldRef()[patchi][faceI] = QcFlame * (mlr/10.0/(Foam::exp(mlr/scalar(10.0))-scalar(1)));
	//	    }
	//	}
	//    }
	//}
 
    }
    
}

// Rarely needed: flame heights written to file...
if (writeFlameHeights && mag(g).value() > 0.0) 
{    
    flameHeight = dimensionedScalar("zero",dimLength,0.0);
    flameHeight2 = dimensionedScalar("zero",dimLength,0.0);
    
    forAll (ft.internalField(),cellI)
    {
	if (ft[cellI] >= ftSt)
	{ 
	    flameHeight[cellI] = mag( cellCentres[cellI] & (g/mag(g)).value() );
	}
	if (ft[cellI] >= 0.2*ftSt)
	{ 
	    flameHeight2[cellI] = mag( cellCentres[cellI] & (g/mag(g)).value() );
	}
    }
    
    outFlameHeight <<  runTime.value() << "  "     
		   << max(flameHeight).value() << "  "
		   << max(flameHeight2).value() << "  "
		   << endl;
}

