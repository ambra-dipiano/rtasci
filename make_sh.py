import sys
import os

script = str(sys.argv[1])
N = int(sys.argv[2])
print(script, N)

rootpath = str(os.path.dirname(os.path.abspath(__file__))).replace('cta-sag-sci', '')
path = f'{rootpath}/jobs'
shname = f'{path}/{script.replace(".py", ".sh")}'
shrun = f'{path}/{script.replace(".py", "_runjobs.sh")}'
sh = open(shname, 'w+')
sh.write('#!/bin/bash\n')
sh.write('\nsource activate scitools')

for i in range(N):
    if i == 0:
        sh.write(f'\n\tpython3 ../{script} 1000 {True}')
        sh.write(f'\n\tpython3 ../{script} 100 {False}')
        sh.write(f'\n\tpython3 ../{script} 10 {False}')
    else:
        sh.write(f'\n\tpython3 ../{script} 1000 {False}')
        sh.write(f'\n\tpython3 ../{script} 100 {False}')
        sh.write(f'\n\tpython3 ../{script} 10 {False}')

sh.write('\n')
sh.close()

sh = open(shrun, 'w+')
sh.write('#!/bin/bash\n')
sh.write(f'\n#SBATCH --job-name={script.replace(".py", "")}')
sh.write(f'\n#SBATCH --output=slurm-$(jobid)-{script.replace(".py", ".out")}')
sh.write(f'\n#SBATCH --account=ambra')
sh.write(f'\n#SBATCH --nodes=1')
sh.write(f'\n#SBATCH --cpus-per-task=1')
sh.write(f'\n\nexec sh ./{script.replace(".py", ".sh")}\n')
sh.close()