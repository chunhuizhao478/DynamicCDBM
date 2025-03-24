Use a single layer of block as damageable region to mimick a fault
Created By Chunhui Zhao, Mar 22th, 2025

#!/bin/bash
#PBS -N buried3d
#PBS -l select=5:ncpus=48:mpiprocs=30
#PBS -l walltime=48:00:00
#PBS -P moose

cd $PBS_O_WORKDIR

module load use.moose moose-dev-openmpi/2025.02.18
mpiexec -n 150 moose-dev-exec ./dynamic_cdbm-opt -i examples/buried_fault/static_solve/static_solve.i
mpiexec -n 150 moose-dev-exec ./dynamic_cdbm-opt -i examples/buried_fault/dynamic_solve/dynamic_solve.i