Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

#include "readAdditionalThermo.H" // kvm

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

#include "createMRF.H"


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

thermo.correct(); // kvm
rho = thermo.rho(); // kvm

mesh.setFluxRequired(p_rgh.name());

#include "phrghEqn.H"


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

volScalarField Qdot_rad
(
    IOobject
    (
        "Qdot_rad",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

volScalarField radT4fac
(
    IOobject
    (
        "radT4fac",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

#include "createDpdt.H"

#include "createK.H"


volScalarField kappa // kvm
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()*thermo.alpha()
);

volScalarField cp // kvm
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()
);

singleStepReactingMixture<gasHThermoPhysics>& singleMixture // kvm
(
    dynamic_cast<singleStepReactingMixture<gasHThermoPhysics>&>
    (thermo)
);

// stoichiometric O2 to fuel ratio
scalar s(singleMixture.s().value()); // kvm

// stoichiometric air to fuel ratio
scalar stoicRatio(singleMixture.stoicRatio().value());   // kvm

// heat of combustion [J/kg]
scalar qFuel(singleMixture.qFuel().value());   // kvm

label fuelIndex(singleMixture.fuelIndex()); // kvm
const label inertIndex(composition.species()[inertSpecie]);
// label inertIndex(singleMixture.inertIndex()); // kvm

const volScalarField* O2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
if (thermo.composition().contains("O2"))  // kvm
{
    O2Ptr = &thermo.composition().Y("O2"); // kvm
}
const volScalarField& O2 = *O2Ptr; // kvm

const volScalarField* CO2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
if (thermo.composition().contains("CO2"))  // kvm
{
    CO2Ptr = &thermo.composition().Y("CO2"); // kvm
}
const volScalarField& CO2 = *CO2Ptr; // kvm

const volScalarField& fu = Y[fuelIndex];  // kvm

scalar YO2Inf = 0.23301; //hardcode for now

// Fuel molecular weight
scalar Wu_ = singleMixture.speciesData()[fuelIndex].W(); // kvm
// Fuel enthalpy of formation
scalar Hu_ = singleMixture.speciesData()[fuelIndex].hc(); // kvm

// compute stoichiometric mixture fraction
scalar ftSt = 1.0 / ( 1.0 + stoicRatio ); // kvm
Info << "stoichiometric mixture fraction is = " << ftSt << nl; // kvm


// create fileds for mixture fraction
volScalarField ft // kvm
(
    IOobject
    (
        "ft",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (fu*s-O2+YO2Inf)/(s+YO2Inf)
);

//- For sprinkler skipping (N.R.)
    volScalarField sprayFluxAF
    (
        IOobject
        (
            "sprayFluxAF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimLength/dimLength/dimTime,0.0)
    );

    volVectorField sprayFluxNet
    (
        IOobject
        (
            "sprayFluxNet",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimMass/dimLength/dimLength/dimTime, Foam::vector(0,0,0)),
        zeroGradientFvPatchVectorField::typeName
    );

    volScalarField sprayTemperature
    (
        IOobject
        (
            "sprayTemperature",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature,273.15),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField sprayMass
    (
        IOobject
        (
            "sprayMass",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass,0.0),
        zeroGradientFvPatchScalarField::typeName
    );


#include "createClouds.H"
#include "createSurfaceFilmModel.H"
#include "createPyrolysisModel.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"
