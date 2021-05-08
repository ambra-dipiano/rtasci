# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import yaml
import sys
import os
import pandas as pd
import argparse

# ---------------------------------------------------------------------------- !

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, default='cfg/myconfig.yml', help='yaml configuration file')
parser.add_argument('--tt', type=float, default=100, help='total trials')
parser.add_argument('--tn', type=float, default=5, help='trials per node')
parser.add_argument('--delay', type=float, default=90, help='delay')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--flux', type=float, default=1, help='flux scaling factor')
parser.add_argument('--env', type=str, default='scitools', help='environment to activate')
parser.add_argument('--pipe', type=str, default='pipe', help='pipeline to run')
args = parser.parse_args()

#print(args)

# compose file path
filename = os.path.join(os.path.expandvars('$PWD'), args.infile)
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
    print(f"run {i+1:02d}")
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
        if args.pipe.lower() == 'pipe':
            f.write(f'\n\tpython rtapipe.py -f {yml}\n')
        elif args.pipe.lower() == 'wilks':
            f.write(f'\n\tpython emptyfields.py -f {yml}\n')


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