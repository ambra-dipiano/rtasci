import yaml
import sys
import os
from os.path import join

# configure
cfgfile = sys.argv[1]
pypath = str(os.path.dirname(os.path.abspath(__file__)))  
configuration = open(os.path.join(pypath, cfgfile))
cfg = yaml.load(configuration, Loader=yaml.FullLoader)

if cfg['options']['extract_data']:
    print('Preparing the catalog...\n')
    os.system(f'python3 simGRBpreparation.py {cfgfile}')

print('\n\nExit\n\n')