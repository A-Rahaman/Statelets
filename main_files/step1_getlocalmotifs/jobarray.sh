#!/bin/bash
#SBATCH -N 1
#SBATCH -n  1
#SBATCH -p  qTRD
#SBATCH -t  7200
#SBATCH -J  mrahaman1
#SBATCH -e error%A.err
#SBATCH -o out%A.out
#SBATCH -A  PSYC0002
#SBATCH --oversubscribe
#SBATCH -c 4
export OMP_NUM_THREADS=1
export MODULEPATH=/apps/Compilers/modules-3.2.10/Debug-Build/Modules/3.2.10/modulefiles/
NODE=$(hostname)
module load Framework/Matlab2016b
matlab -r 'getshapeletsPairWise($SLURM_ARRAY_TASK_ID)' -nodisplay
