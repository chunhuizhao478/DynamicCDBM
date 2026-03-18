#!/bin/bash
#SBATCH -J zip_on_fault                # Job name
#SBATCH -o zip_on_fault.o%j            # Name of stdout output file
#SBATCH -e zip_on_fault.e%j            # Name of stderr error file
#SBATCH -p small                       # Queue (partition) name
#SBATCH -N 1                           # Total # of nodes
#SBATCH -n 56                          # Total # of tasks
#SBATCH -t 01:00:00                    # Run time (hh:mm:ss)
#SBATCH --mail-type=all                # Send email at begin and end of job
#SBATCH -A EAR20006                    # Project/Allocation name
#SBATCH --mail-user=chunhui3@illinois.edu

BASEDIR=/scratch2/10024/zhaochun/projects/DynamicCDBM_03162026/examples/3dslipweakening/parametric_study/dynamic_code_files

cd $BASEDIR

tar cf - dynamic_solve_alpha0d3_dsigma13_Dc0d4_db_csv_main_fault_* | pigz -p 56 > alpha0d3_dsigma13_Dc0d4_db_on_fault.tar.gz

echo "Zipping complete."
