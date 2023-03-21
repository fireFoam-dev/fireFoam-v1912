#!/bin/bash

# Process options
VERSION=1912
OPENFOAM_BRANCH=maintenance-v$VERSION
BUILD_TYPE=Opt # Opt or Debug
START_STEP=1

# External links
THIRDPARTY_LINK=https://sourceforge.net/projects/openfoam/files/v$VERSION/ThirdParty-v$VERSION.tgz
OPENFOAM_LINK=https://develop.openfoam.com/Development/openfoam.git
FIREFOAM_LINK=https://github.com/fireFoam-dev/fireFoam-v$VERSION.git

# ----- No changes below this line ------ #

# Set installation directory and local (cloned) repository names
MYINSTDIR=$PWD
MYFIREFOAM_REPO=fireFoam-v$VERSION
MYOPENFOAM_REPO=OpenFOAM-v$VERSION

function setup_environment() {
  echo ========================================
  echo "Setting up environment"
  source $MYINSTDIR/$MYOPENFOAM_REPO/etc/bashrc WM_COMPILE_OPTION=$BUILD_TYPE
  echo "...done."
  echo ========================================
}

echo "Beginning build on $(date)"

#1) Get ThirdParty package 
if [ $START_STEP -eq 1 ]
then
  cd $MYINSTDIR
  echo ========================================
  echo "Getting ThirdParty package"
  wget $THIRDPARTY_LINK
  if [ ! -f ThirdParty-v$VERSION.tgz ]
  then
    echo "> Downloading ThirdParty failed!"
    exit 2
  fi
  tar -zxvf ThirdParty-v$VERSION.tgz
  rm ThirdParty-v$VERSION.tgz
  echo "...done."
  echo ========================================
fi

#2) Clone OpenFOAM repository
if [ $START_STEP -le 2 ]
then
  cd $MYINSTDIR
  echo ========================================
  echo "Cloning OpenFOAM repository"
  git clone $OPENFOAM_LINK $MYOPENFOAM_REPO
  cd $MYOPENFOAM_REPO
  git submodule init
  echo "...done."
  echo ========================================
fi

#3) Clone FireFOAM repository
if [ $START_STEP -le 3 ]
then
  cd $MYINSTDIR
  echo ========================================
  echo "Cloning FireFOAM repository"
  git clone $FIREFOAM_LINK $MYFIREFOAM_REPO
  echo "...done."
  echo ========================================
fi

#4) Checkout OpenFOAM branch needed for FireFOAM
if [ $START_STEP -le 4 ]
then
  cd $MYINSTDIR/$MYOPENFOAM_REPO
  echo ========================================
  echo "Checking out OpenFOAM branch: $OPENFOAM_BRANCH"
  git checkout $OPENFOAM_BRANCH
  echo "...done."
  echo ========================================
fi

#5) Compile Third Party software
if [ $START_STEP -le 5 ]
then
  setup_environment
  cd $WM_THIRD_PARTY_DIR
  echo ========================================
  echo "Compiling ThirdParty..."
  ./Allwmake -j -l
  echo "...done."
  echo ========================================
fi

#6) Compile OpenFOAM
if [ $START_STEP -le 6 ]
then
  setup_environment
  cd $WM_PROJECT_DIR
  echo ========================================
  echo "Compiling OpenFOAM..."
  ./Allwmake -j -l
  echo "...done."
  echo ========================================
fi

#7) Build fireFoam
if [ $START_STEP -le 7 ]
then
  setup_environment
  cd $MYINSTDIR/$MYFIREFOAM_REPO
  echo ========================================
  echo "Compiling FireFOAM..."
  ./Allwmake -j -l
  echo "...done."
  echo ========================================
fi

echo "Build completed on $(date)"
