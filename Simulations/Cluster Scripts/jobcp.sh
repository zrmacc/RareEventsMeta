#!/bin/bash
#SBATCH --account=def-jgrons
#SBATCH --array=1-40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=0-06:00
#SBATCH --job-name=meta24
#SBATCH --output=meta24-%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ysyyjwh@gmail.com

module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.1
module load openmpi/3.1.4

echo "MY slurm arrary task id: " $SLURM_ARRAY_TASK_ID

echo "MY slurm array job id: " $SLURM_ARRAY_JOB_ID

Rscript Rscripts/cp_array_NA.R --studies 24 --rate 0.06 --alpha 140 --beta 70 --reps 50 --mc 2000 --bl 1 --job ${SLURM_ARRAY_TASK_ID} --out "Results/K24/"



