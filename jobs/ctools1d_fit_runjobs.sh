
#SBATCH --job-name=ctools1d_fit
#SBATCH --output=slurm-$(jobid)-ctools1d_fit.out
#SBATCH --account=ambra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

exec sh /home/ambra/Desktop/CTA/projects/cta-sag-sci/jobs/ctools1d_fit.sh
