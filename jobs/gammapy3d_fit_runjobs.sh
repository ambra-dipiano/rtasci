
#SBATCH --job-name=gammapy3d_fit
#SBATCH --output=slurm-$(jobid)-gammapy3d_fit.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh /home/ambra/Desktop/CTA/projects/cta-sag-sci/jobs/gammapy3d_fit.sh
