#!/bin/bash

#SBATCH --job-name=gammapy3d_fit
#SBATCH --output=slurm-$(jobid)-gammapy3d_fit.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh ./gammapy3d_fit.sh
