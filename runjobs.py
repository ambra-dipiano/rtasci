from os import listdir, system
from os.path import isfile, join

path = '/home/ambra/Desktop/CTA/projects/cta-sag-sci/jobs/'

jobs = [f for f in listdir(path) if '_runjobs.sh' in f and isfile(join(path, f))]

for job in jobs:
    #os.system('sbatch job')
    print(f'sbatch {job}')
