# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import argparse
from rtasci.cfg.Config import Config

# configure
parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
parser.add_argument('--merge', type=str, default='true', help='Merge in single phlist (true) or use observation library (false)')
parser.add_argument('--remove', type=str, default='true', help='Keep only outputs')
parser.add_argument('--print', type=str, default='false', help='Print out results')
args = parser.parse_args()

cfg = Config(args.cfgfile)

# simulations
if cfg.get('simtype').lower() == 'grb':
    if cfg.get('extract_data'):
        print('\nPreparing GRB catalog...\n')
        os.system(f'python3 prepareGRBcatalog.py -f {args.cfgfile}')
    print('\nRun simulations...\n')
    os.system(f'python3 simGRBcatalog.py -f {args.cfgfile} --merge {args.merge.lower()} --remove {args.remove.lower()} --print {args.print.lower()}')
elif cfg.get('simtype').lower() == 'bkg':
    print('\nComputing BKG-ONLY simulations is work in progress\n')
elif cfg.get('simtype').lower() == 'wilks':
    print("\nRun empty fields simulation + analysis")
    os.system(f"python3 emptyfields.py -f {args.cfgfile} --remove {args.remove.lower()} --print {args.print.lower()}")
elif cfg.get('simtype').lower() == 'skip':
    pass
else:
    raise ValueError('Invalid "simtype" value')

# analysis
if cfg.get('tool') not in ('ctools', 'gammapy', 'rtatool'):
    raise ValueError('Invalit "tool" selection.')
else:
    print(f'\nRun analysis...\n')
    pipeline = f"{cfg.get('tool')}{cfg.get('type')}"
    if cfg.get('blind'):
        pipeline += '_blind'
    if not cfg.get('binned'):
        pipeline += '_unbinned'
    pipeline += '.py'
    print(f'Pipeline: {pipeline}')
    os.system(f"python3 pipelines/{pipeline} -f {args.cfgfile} --merge {args.merge.lower()} --remove {args.remove.lower()} --print {args.print.lower()}")

if "_trials" in args.cfgfile:
    os.system(f"rm {args.cfgfile}")
print('\n\nExit\n\n')