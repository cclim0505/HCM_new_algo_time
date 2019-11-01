#!/bin/bash
#PBS -N haus
#PBS -q serial
#PBS -l select=1:ncpus=1
#PBS -l walltime=96:00:00
#PBS -P MST108285
#PBS -o PBS.log
#PBS -e PBS.err
####PBS -r n
cd $PBS_O_WORKDIR
module load intel/2018_u1
./Hausdorff.out
