    dictionary dictMULES =
    mesh.solutionDict().subDict("MULES");

    Switch isMULES(dictMULES.lookup("bMULES"));
    Switch isHsMULES(dictMULES.lookupOrDefault("bHsMULES", false));

