#!/bin/sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions
# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase

runApplication foamCleanPolyMesh
mv log.foamCleanPolyMesh log.foamCleanPolyMesh.gasPhase
runApplication foamCleanPolyMesh -region filmRegion
mv log.foamCleanPolyMesh log.foamCleanPolyMesh.filmRegion
runApplication foamCleanPolyMesh -region cardboardBoxRegion
mv log.foamCleanPolyMesh log.foamCleanPolyMesh.cardboardBox
runApplication foamCleanPolyMesh -region woodPalletRegion
mv log.foamCleanPolyMesh log.foamCleanPolyMesh.woodPallet

rm -fr log*
rm -fr *.obj
rm -rf VTK
rm -rf patch*
rm -rf fieldMimMax  HRR  probes
rm -fr fieldMinMax
rm -fr postProcessing

# Run
#runParallel $application 4

# -----------------------------------------------------------------------------
