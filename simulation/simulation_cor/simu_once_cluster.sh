#!/bin/bash

#SBATCH --array=1-1000
#SBATCH --nodes=1
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1

## turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx

ml R/4.4.2-gfbf-2024a 
R CMD BATCH --no-save simu_once_cluster.R Rout/simu_once_cluster_${SLURM_ARRAY_TASK_ID}.Rout