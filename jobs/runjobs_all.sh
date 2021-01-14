#!/bin/bash

#SBATCH --job-name=allscripts
#SBATCH --output=slurm-allscripts.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/gammapy1d_binned_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/gammapy3d_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/ctools1d_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/gammapy3d_binned_blindfit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/ctools3d_binned_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/ctools3d_blindfit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/slurm-ctools1d_binned_fit.out
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/ctools3d_binned_blindfit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/gammapy3d_blindfit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/ctools3d_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/gammapy1d_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/gammapy3d_binned_fit.sh
exec sh /data01/homes/cta/gammapy_integration/cta-sag-sci/jobs/ctools1d_binned_fit.sh
