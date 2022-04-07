IOdictionary additionalControlsDict
(
    IOobject
    (
        "additionalControls",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

bool solvePrimaryRegion
(
    additionalControlsDict.lookupOrDefault("solvePrimaryRegion", true)
);

bool solvePyrolysisRegion
(
    additionalControlsDict.lookupOrDefault("solvePyrolysisRegion", true)
);

scalar maxDi = pyrolysis.maxDiff();

// Compute and write outFlameHeights_<startTime> file
const bool writeFlameHeights
(
    additionalControlsDict.lookupOrDefault("writeFlameHeights", false)
);
// Output frequency of infoOutput
const word infoOutputControl
(
    additionalControlsDict.lookupOrDefault<word>("infoOutputControl", "none")
);
const label infoOutputInterval
(
    additionalControlsDict.lookupOrDefault<label>("infoOutputInterval", 1)
);
// Compute diagnostics
const bool writeInfoOutput
(
    infoOutputControl == "timeStep"
);

// Update mesh - needed for dynamicFvMesh 
const bool updateMesh
(
    additionalControlsDict.lookupOrDefault("updateMesh",false)
);

const label updateMeshFrequency
(
    additionalControlsDict.lookupOrDefault("updateMeshFrequency",-1)
);


// Also compute the rarely used diagnostics
const bool verboseInfoOutput
(
    writeInfoOutput
 && additionalControlsDict.lookupOrDefault("verboseInfoOutput", false)
);
// Skip the computations actually used in the solver (Use carefully: Could lead to solution errors)
const bool skipInfoOutputSolverCalcs
(
    additionalControlsDict.lookupOrDefault("skipInfoOutputSolverCalcs", false)
);
if ( skipInfoOutputSolverCalcs )
{
    WarningIn
    (
	"createRegionControls.H"
    )   << "\nYOU HAVE DISABLED THE CALCULATION OF THE FOLLOWING FIELDS: \n\n" 
	<< "\t    convectiveHeatFlux_T \n"
	<< "\t    rhoSpray \n"
	<< "\t    mlrSolid \n\n"
	<< "THIS WILL CAUSE SOLUTION ERRORS IF THESE VARIABLES ARE NEEDED IN YOUR SIMULATIONS!\n\n"
	<< "TO ENABLE THE COMPUTATIONS, SET 'skipInfoOutputSolverCalcs false' IN 'constant/additionalControls'\n"
	<< "...OR SIMPLY DELETE THE KEYWORD ENTRY\n\n"
	<< endl;
}
