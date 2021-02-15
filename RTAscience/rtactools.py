import yaml
import sys
import os
import argparse
from os.path import join
from RTAscience.cfg.Config import Config

# configure
parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
args = parser.parse_args()

cfg = Config(args.cfgfile)

if cfg.get('simtype').lower() == 'grb':
    if cfg.get('extract_data'):
        print('\nPreparing GRB catalog...\n')
        os.system(f'python3 simGRBpreparation.py -f {args.cfgfile}')
    print('\nRun simulations...\n')
    os.system(f'python3 simGRBcat.py -f {args.cfgfile}')
elif cfg.get('simtype').lower() == 'bkg':
    print('\nComputing BKG-ONLY simulations is work in progress\n')
else:
    raise ValueError('Invalid "simtype" value')

print('\n\nExit\n\n')