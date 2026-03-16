#!/bin/bash
#SBATCH -J dynamic_solve_alpha0_dsigma7_Dc0d4_db        # Job name
#SBATCH -o dynamic_solve_alpha0_dsigma7_Dc0d4_db.o%j    # Name of stdout output file
#SBATCH -e dynamic_solve_alpha0_dsigma7_Dc0d4_db.e%j    # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 30              # Total # of nodes
#SBATCH -n 400             # Total # of mpi tasks
#SBATCH -t 30:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A EAR20006        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=chunhui3@illinois.edu

# Load necessary modules
module swap intel gcc
#module swap impi mvapich2-x
module load cuda
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Set compilers and flags
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77

export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Enable MPI debugging
export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

export MOOSE_JOBS=6 METHODS=opt

ibrun /scratch1/10024/zhaochun/projects/DynamicCDBM_03162026/dynamic_cdbm-opt -i /scratch1/10024/zhaochun/projects/DynamicCDBM_03162026/examples/3dslipweakening/parametric_study/dynamic_code_files/dynamic_solve_alpha0_dsigma7_Dc0d4_db.i --allow-unused
