#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication topoSet
runApplication createPatch -overwrite

#set orientation of horizontal sampling planes
nPlanes=29
for ((i=1; i<=9; ++i))
do
        j=$( echo "scale=1; $i/10" | bc )
        orientFaceZone "plane_z0$j" "(0 0 10)"
done

for ((i=10; i<=$nPlanes; ++i))
do
        j=$( echo "scale=1; $i/10" | bc )
        orientFaceZone "plane_z$j" "(0 0 10)"
done


cp 0/ph_rgh-orig 0/ph_rgh

runApplication decomposePar -force

#------------------------------------------------------------------------------
