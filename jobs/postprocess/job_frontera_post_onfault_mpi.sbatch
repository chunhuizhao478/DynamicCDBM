#!/bin/bash
#SBATCH -J post_onfault_mpi            # Job name
#SBATCH -o post_onfault_mpi.o%j        # Name of stdout output file
#SBATCH -e post_onfault_mpi.e%j        # Name of stderr error file
#SBATCH -p development                 # Queue (partition) name
#SBATCH -N 1                           # Total # of nodes
#SBATCH -n 56                          # Total # of tasks
#SBATCH -t 02:00:00                    # Run time (hh:mm:ss)
#SBATCH --mail-type=all                # Send email at begin and end of job
#SBATCH -A EAR20006                    # Project/Allocation name
#SBATCH --mail-user=chunhui3@illinois.edu

module load python3
pip install --user tqdm

SCRIPT=/scratch2/10024/zhaochun/projects/DynamicCDBM_03162026/examples/3dslipweakening/parametric_study/postprocess/post_onfault_mpi.py

ibrun python3 $SCRIPT

echo "Post-processing complete."
