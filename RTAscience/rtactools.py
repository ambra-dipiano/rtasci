import yaml
import sys
import os
import argparse
from os.path import join
from RTAscience.cfg.Config import Config

# configure
parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
parser.add_argument('--remove', type=str, default='true', help='Keep only outputs')
args = parser.parse_args()

cfg = Config(args.cfgfile)

# simulations
if cfg.get('simtype').lower() == 'grb':
    if cfg.get('extract_data'):
        print('\nPreparing GRB catalog...\n')
        os.system(f'python3 prepareGRBcatalog.py -f {args.cfgfile}')
    print('\nRun simulations...\n')
    os.system(f'python3 simGRBcatalog.py -f {args.cfgfile}')
elif cfg.get('simtype').lower() == 'bkg':
    print('\nComputing BKG-ONLY simulations is work in progress\n')
elif cfg.get('simtype').lower() == 'skip':
    pass
else:
    raise ValueError('Invalid "simtype" value')

# analysis
if cfg.get('tool') not in ('ctools', 'gammapy', 'rtatool'):
    raise ValueError('Invalit "tool" selection.')
elif cfg.get('tool') != 'ctools':
    raise ValueError('Option not yet implemented.')
elif cfg.get('tool') == 'ctools':
    print(f'\nRun analysis...\n')
    pipeline = f"{cfg.get('tool')}{cfg.get('type')}_"
    if cfg.get('binned'):
        pipeline += 'binned_'
    if cfg.get('blind'):
        pipeline += 'blind'
    pipeline += 'fit.py'
    print(f'Pipeline: {pipeline}')
    os.system(f"python3 pipelines/{pipeline} -f {args.cfgfile} --remove {args.remove.lower()}")

print('\n\nExit\n\n')