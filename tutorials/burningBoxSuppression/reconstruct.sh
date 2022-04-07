#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

runApplication -s primary reconstructPar 
runApplication -s film    reconstructPar -region filmRegion 
runApplication -s fuel    reconstructPar -region fuelRegion 


exit;
