import yaml
import sys
import os
from os.path import join

# configure
cfgfile = sys.argv[1]

print('Preparing the catalog...\n')
os.system(f'python3 simGRBpreparation.py {cfgfile}')

print('\n\nExit\n\n')