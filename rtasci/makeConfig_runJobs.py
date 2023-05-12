# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
# *******************************************************************************

import os
import yaml
import argparse
from pathlib import Path

# ---------------------------------------------------------------------------- !

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, required=True, help='yaml configuration file')
parser.add_argument('--tt', type=int, required=True, help='total trials')
parser.add_argument('--cpus', type=int, required=True, help='number of cpus')
parser.add_argument('--script', type=str, required=True, help='script to run')
parser.add_argument('--env', type=str, required=True, help='environment to activate')
parser.add_argument('--delay', type=float, default=90, help='delay')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--flux', type=float, default=1, help='flux scaling factor')
parser.add_argument('--print', type=str, default='false', help='print checks and outputs')
parser.add_argument('-out', '--output-dir', type=str, required=False, default="", help='The path to the output directory')
args = parser.parse_args()

if "DATA" not in os.environ:
    raise ValueError("Please, export DATA")
print(f'\nDATA={os.environ["DATA"]}')

# compose file path
filename = Path(args.infile)
if not filename.is_file():
    raise ValueError('configuration file not found')
# load yaml
with open(filename) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

output_path = filename.parent
ext = filename.suffix
name = filename.stem

trials_per_cpu = int(args.tt / args.cpus)

if args.off == 'gw':
    config['simulation']['offset'] = args.off
else:
    config['simulation']['offset'] = float(args.off)
config['simulation']['delay'] = args.delay
config['setup']['trials'] = int(trials_per_cpu)
config['setup']['scalefluxfactor'] = args.flux
config['options']['plotsky'] = False
start_count = config['setup']['start_count']

print(f"SLURM configuration:\n\tNumber of jobs: {args.cpus}\n\tTotal trials: {args.cpus*trials_per_cpu}\n\tTrials per job: {trials_per_cpu}\n\tStart count: {start_count}\n\tOutput dir: {args.output_dir}")
input("Press any key to start!")

for i in range(args.cpus):
    run = f"job_{i+1:02d}"
    job_name = f'trials_{i*trials_per_cpu+1}-{(i+1)*trials_per_cpu}'
    #print("\n")
    #print(f"run: {run}")
    #print(f"job name: {job_name}")
    output_dir = output_path.joinpath(run)
    output_dir.mkdir(parents=True, exist_ok=True)

    # save new config
    config['setup']['start_count'] = start_count + i*trials_per_cpu
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

        script_name = Path(args.script).stem.lower()
        if script_name == 'rtapipe':
            f.write(f'\n\tpython {args.script}.py -f {config_outname} --print {args.print.lower()}\n')
        elif script_name == 'wilks':
            f.write(f'\n\tpython {args.script} -f {config_outname}\n')
        elif script_name == "simbkg":
            f.write(f'\n\tpython {args.script} -f {config_outname} -out {args.output_dir}\n')
        else:
            raise ValueError(f"Script {script_name} is not supported.")

    # write job
    job_outname = output_dir.joinpath(f"job_{job_name}").with_suffix(".sh")
    job_outlog = output_dir.joinpath(f"job_{job_name}").with_suffix(".log")

    with open(job_outname, 'w+') as f:
        f.write('#!/bin/bash')
        f.write(f'\n\n#SBATCH --job-name={script_name}-slurm-job')
        f.write(f'\n#SBATCH --output={job_outlog}')
        f.write('\n#SBATCH --account=baroncelli')
        f.write('\n#SBATCH --partition=large_lc')
        f.write(f'\n\nexec sh {str(sh_outname)}\n')

    #print(f"Configuration file={config_outname}")
    #print(f"Exec file created={sh_outname}")
    #print(f"Slurm configuration file created={job_outname}")

    os.system(f'sbatch {job_outname}')
