#!/bin/bash

#SBATCH --job-name=gammapy1d_fit
#SBATCH --output=slurm-$(jobid)-gammapy1d_fit.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh ./gammapy1d_fit.sh
