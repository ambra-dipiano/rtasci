#!/bin/bash

#SBATCH --job-name=ambra
#SBATCH --output=slurm-time_gammapy3d_blindfit.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/RTAscience/jobs/time_gammapy3d_blindfit.sh
