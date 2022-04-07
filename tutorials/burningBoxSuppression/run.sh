#!/bin/bash
#SBATCH --job-name=UUP-block
#SBATCH --partition=nathaz_haswell
#SBATCH --threads-per-core=1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=./run-%j
#SBATCH --no-requeue

mpirun -np 8 fireFoam -parallel
