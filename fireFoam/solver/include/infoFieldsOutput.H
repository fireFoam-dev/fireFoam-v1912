    // ft-place-holder

    // create fileds for patch integration of HRR
    surfaceScalarField HRR_fu
    ( 
        IOobject
        (
            "HRR_fu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(fu)*qFuel
    );

    // for outputing flame height
    OFstream outFlameHeight("outFlameHeight_"+runTime.timeName());
    outFlameHeight << "# Time    " << "FlameHeight  "  << endl;

    volScalarField flameHeight
    (
        IOobject
        (
            "flameHeight",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );
    volScalarField flameHeight2
    (
        IOobject
        (
            "flameHeight2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );

    const pointField& cellCentres = mesh.cellCentres();


    // for check ft conservation
    surfaceScalarField phiFt
    ( 
        IOobject
        (
            "phiFt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    surfaceScalarField phiO2
    ( 
        IOobject
        (
            "phiO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    surfaceScalarField phiCO2
    ( 
        IOobject
        (
            "phiCO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    // create fileds for surface integration of sensible enthalpy flux
    surfaceScalarField phiHs
    (
        IOobject
        (
            "phiHs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(thermo.he()) 
    );

    // create fileds for surface integration of chemical enthalpy flux
    volScalarField hc(thermo.hc()); 
    surfaceScalarField phiHc
    (
        IOobject
        (
           "phiHc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create fileds for surface integration of total enthalpy flux
    volScalarField h(hc + thermo.he());
    surfaceScalarField phiH
    (
        IOobject
        (
            "phiH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create UT field for average purpose
    volVectorField UT
    (
        IOobject
        (
            "UT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U*T
    );

    // create rhoU field for average purpose
    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    surfaceScalarField momentumX
    (
        IOobject
        (
            "momentumX",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("momentumX", dimMass*dimLength/pow(dimTime,2), scalar(0))
    );

    surfaceScalarField momentumY
    (
        IOobject
        (
            "momentumY",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("momentumY", dimMass*dimLength/pow(dimTime,2), scalar(0))
    );

    surfaceScalarField momentumZ
    (
        IOobject
        (
            "momentumZ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("momentumZ", dimMass*dimLength/pow(dimTime,2), scalar(0)) // momentum is kg m /s
    );
        //dimensionedScalar("momentumZ", dimMass*dimLength/dimTime, scalar(0))
        //dimensionedScalar("momentumZ", dimMass/(dimLength*dimTime), scalar(0))
 //dimensionedScalar unitMom
 //(
 //    "unitMom",
 //    dimTime,
 //    1.0
 //);


    // create rhoT field for average purpose
    volScalarField rhoT
    (
        IOobject
        (
            "rhoT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*T
    );

    // create rhoFu field for average purpose
    volScalarField rhoFu
    (
        IOobject
        (
            "rhoFu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*fu
    );

    // create B field for average purpose
    volSymmTensorField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(U)
    );


    surfaceScalarField convectiveHeatFlux_L
    (
        IOobject
        (
            "convectiveHeatFlux_L",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

    surfaceScalarField convectiveHeatFlux_T
    (
        IOobject
        (
            "convectiveHeatFlux_T",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

/*
    volScalarField wallConvectiveHeatFlux
    (
        IOobject
        (
            "wallConvectiveHeatFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );
*/

    volScalarField QcWallFunction
    (
        IOobject
        (
            "QcWallFunction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

    IOdictionary QcWallFunctionDict
    (
        IOobject
        (
            "QcWallFunctionDict",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            //IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar QcFlame =
    QcWallFunctionDict.lookupOrDefault<scalar>("QcFlame", 16000.0);
    
    scalar QcThreshold =
    QcWallFunctionDict.lookupOrDefault<scalar>("QcThreshold", 5000.0);

    volScalarField alphaM
    (
        IOobject
        (
            "alphaM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->alpha()
    );

    volScalarField muM
    (
        IOobject
        (
            "muM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->mu()
    );

    volScalarField d
    (
        IOobject
        (
            "d",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->alpha()/lewisNo/rho
    );

    volScalarField dSgs
    (
        IOobject
        (
            "dSgs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (turbulence->alphaEff()-turbulence->alpha())/rho
    );

    //- Flame Extinction Model (Bert)
    volScalarField volFracSpray
    (
        IOobject
	(
            "VFSpray",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    volScalarField rhoSpray
    (
        IOobject
	(
            "rhoSpray",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", rho.dimensions(),1000.0)
    );

    dimensionedScalar rhoWater("rhoWater",rho.dimensions(),1000.0);

    volScalarField mlrSolid
    (
        IOobject
	(
            "mlrSolid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimArea/dimTime,0.0)
    );
