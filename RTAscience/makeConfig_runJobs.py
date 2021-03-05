# RUNNING THE SCRIPT:
# $ python configYamlFilters.py CONFIGURATION.yaml YEARS_args.indir
# example:
# $ python configYamlFilters.py M10B_294278402_296870402.yaml YEARS5g019

import yaml
import sys
import os
import pandas as pd
import argparse

# ---------------------------------------------------------------------------- !

parser = argparse.ArgumentParser()
parser.add_argument('--infile', type=str, default='cfg/myconfig.yml', help='yaml configuration file')
parser.add_argument('--tt', type=float, default=10000, help='total trials')
parser.add_argument('--tn', type=float, default=500, help='trials per node')
parser.add_argument('--delay', type=float, default=50, help='delay')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--flux', type=float, default=1, help='flux scaling factor')
parger.add_argument('--env', type=str, default='ctools', help='environment to activate')
args = parser.parse_args()

#print(args)

# compose file path
filename = os.path.join(os.path.expandvars('$PWD'), args.infile)
print(filename)
if not os.path.isfile(filename):
    raise ValueError('yaml file not found')
# load yaml
with open(filename) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

if args.off == 'gw':
    config['simulation']['offset'] = args.off
else:
    config['simulation']['offset'] = float(args.off)
config['simulation']['delay'] = args.delay
config['setup']['trials'] = int(args.tn)
config['setup']['scalefluxfactor'] = args.flux
config['options']['plotsky'] = False

for i in range(int(args.tt/args.tn)):
    # save new gonfig
    config['setup']['start_count'] = int(i*args.tn)
    outname = args.infile.replace('.yml',f'_trials{i*args.tn+1}-{(i+1)*args.tn}')
    yml = outname + '.yml' 
    if os.path.isfile(yml):
        os.remove(yml)
    with open(yml, 'w+') as f:
        new_config = yaml.dump(config, f, default_flow_style=False)

    # write bash
    sh = outname.replace('cfg/', 'jobs/') + '.sh'
    with open(sh, 'w+') as f:
        f. write('#!/bin/bash\n')
        f.write(f'\nsource activate {args.env}')
        f.write('\n\texport DATA=/data01/homes/cta/gammapy_integration/DATA/')
        f.write(f'\n\tpython rtactools.py -f {yml}\n')

    """ 
    # write job
    job = outname.replace('cfg/', 'jobs/job_') + '.sh'
    with open(job, 'w+') as f:
        f.write('#!/bin/bash')
        f.write('\n\n#SBATCH --job-name=' + outname)
        f.write('\n#SBATCH --output=slurm-' + outname+ '.out')
        f.write('\n#SBATCH --account=rt')
        f.write('\n#SBATCH --ntasks=1')
        f.write('\n#SBATCH --nodes=1')
        f.write('\n#SBATCH --cpus-per-task=1')
        f.write('\n\nexec sh ' + sh + '\n') """

    os.system(f'sbatch {sh}')