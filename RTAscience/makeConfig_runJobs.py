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
from pathlib import Path
import pandas as pd
import argparse

# ---------------------------------------------------------------------------- !

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, required=True, help='yaml configuration file')
parser.add_argument('--tt', type=int, required=True, help='total trials')
parser.add_argument('--tn', type=int, required=True, help='trials per node')
parser.add_argument('--script', type=str, required=True, help='script to run')
parser.add_argument('--env', type=str, required=True, help='environment to activate')
parser.add_argument('--delay', type=float, default=90, help='delay')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--flux', type=float, default=1, help='flux scaling factor')
parser.add_argument('--print', type=str, default='false', help='print checks and outputs')
args = parser.parse_args()

if "DATA" not in os.environ:
    raise ValueError("Please, export DATA")
    exit(0)
print(f'\nDATA={os.environ["DATA"]}')

# compose file path
#filename = os.path.join(os.path.expandvars('$PWD'), args.infile)
filename = Path(args.infile)
if not filename.is_file():
    raise ValueError('configuration file not found')
# load yaml
with open(filename) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

output_path = filename.parent
ext = filename.suffix
name = filename.stem

if args.off == 'gw':
    config['simulation']['offset'] = args.off
else:
    config['simulation']['offset'] = float(args.off)
config['simulation']['delay'] = args.delay
config['setup']['trials'] = int(args.tn)
config['setup']['scalefluxfactor'] = args.flux
config['options']['plotsky'] = False

number_of_jobs=int(args.tt/args.tn)
print(f"SLURM configuration:\n\tNumber of jobs: {number_of_jobs}\n\tTotal trials: {args.tt}\n\tTrials per job: {args.tn}")
input("Press any key to start!")

for i in range(number_of_jobs):
    run = f"job_{i+1:02d}"
    job_name = f'trials_{i*args.tn+1}-{(i+1)*args.tn}'
    #print("\n")
    #print(f"run: {run}")
    #print(f"job name: {job_name}")
    output_dir = output_path.joinpath(run)
    output_dir.mkdir(parents=True, exist_ok=True)

    # save new config
    config['setup']['start_count'] = int(i*args.tn)
    config_outname = output_dir.joinpath(f"{name}_{job_name}").with_suffix(ext)

    if config_outname.is_file():
        config_outname.unlink()
    with open(config_outname, 'w+') as f:
        new_config = yaml.dump(config, f, default_flow_style=False)

    # write bash
    sh_outname = output_dir.joinpath(f"sh_{job_name}").with_suffix(".sh")
    with open(sh_outname, 'w+') as f:
        f. write('#!/bin/bash\n')
        f.write(f'\nsource activate {args.env}')
        f.write(f'\n\texport DATA={os.environ["DATA"]}')

        scriptName = Path(args.script).stem.lower()

        if scriptName == 'pipe':
            f.write(f'\n\tpython {args.script} -f {config_outname} --print {args.print.lower()}\n')
        elif scriptName == 'wilks':
            f.write(f'\n\tpython {args.script} -f {config_outname}\n')
        elif scriptName == "simbkg":
            f.write(f'\n\tpython {args.script} -f {config_outname}\n')
        else:
            raise ValueError(f"Script {scriptName} is not supported.")

    # write job
    job_outname = output_dir.joinpath(f"job_{job_name}").with_suffix(".sh")
    job_outlog = output_dir.joinpath(f"job_{job_name}").with_suffix(".log")

    with open(job_outname, 'w+') as f:
        f.write('#!/bin/bash')
        f.write(f'\n\n#SBATCH --job-name=CTA-sim-slurm-job_{job_name}')
        f.write(f'\n#SBATCH --output={job_outlog}')
        f.write('\n#SBATCH --account=baroncelli')
        f.write('\n#SBATCH --ntasks=1')
        f.write('\n#SBATCH --nodes=1')
        f.write('\n#SBATCH --cpus-per-task=1')
        f.write(f'\n\nexec sh {str(sh_outname)}\n')

    #print(f"Configuration file={config_outname}")
    #print(f"Exec file created={sh_outname}")
    #print(f"Slurm configuration file created={job_outname}")

    os.system(f'sbatch {job_outname}')
