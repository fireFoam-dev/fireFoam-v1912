#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

./clean.sh

runApplication blockMesh

runApplication snappyHexMesh -overwrite

runApplication -s burner topoSet -dict system/topoSetDict.burner

runApplication createPatch -overwrite

runApplication -s pyrolsis topoSet -dict system/topoSetDict.pyrolysis

transformPoints -scale '(0.0254 0.0254 0.0254)'

runApplication -s film extrudeToRegionMesh -overwrite -dict system/extrudeToRegionMeshDict.film

# change boundary patch name
sed -i 's/samplePatch     region0_to_filmRegion_fuel;/samplePatch     region0_to_fuelRegion_fuel;/g' constant/filmRegion/polyMesh/boundary

runApplication -s fuel extrudeToRegionMesh -overwrite -dict system/extrudeToRegionMeshDict.fuel

runApplication -s primary decomposePar -decomposeParDict system/decomposeParDict -force
runApplication -s film    decomposePar -region filmRegion -decomposeParDict system/decomposeParDict
runApplication -s fuel    decomposePar -region fuelRegion -decomposeParDict system/decomposeParDict

cp 0/ph_rgh.orig 0/ph_rgh

exit;
