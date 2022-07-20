#!/bin/bash
#SBATCH --job-name=tutorial
#SBATCH --partition=fire_haswell
#SBATCH --threads-per-core=1
#SBATCH --nodes=8
#SBATCH --ntasks=128
#SBATCH --output=./slurm.%j.out

source $HOME/codes/public_FireFOAM/OpenFOAM-v1912/etc/bashrc
mpirun -np 128 fireFoam -parallel 
