#!/bin/bash

#SBATCH --job-name=ctools3d_blindfit
#SBATCH --output=slurm-$(jobid)-ctools3d_blindfit.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh ./ctools3d_blindfit.sh
