import os
import pandas as pd
import numpy as np
from lib.RTAVisualise import *
from os.path import join, isdir, isfile
from os import listdir

pypath = str(os.path.dirname(os.path.abspath(__file__)))  
datapath = pypath.replace('cta-sag-sci', 'DATA/outputs/crab')
pngpath = join(datapath, 'png')
if not isdir(pngpath):
    print('Creating png folder...')
    os.mkdir(pngpath)

tables = [f for f in listdir(datapath) if 'ctools1d_binned' in f]
print(f'Tables: {len(tables)}\n')

for i, table in enumerate(tables):
    print(f'Collecting data from {table}')
    data = pd.read_csv(join(datapath, table), sep=' ')
    print(data[:5])

    d1000 = data[data['texp'] == 1000]
    d100 = data[data['texp'] == 100]
    d10 = data[data['texp'] == 10]

    m1000 = d1000.mean()
    print(m1000.keys())
    print(m1000[m1000.keys()[0]])

    sci = ['texp', 'sqrt_ts', 'flux', 'flux_err']
    time = ['ttotal', 'timport', 'tsetup', 'tobs', 'tconf', 'tred', 'tblind', 'tcube', 'texpcube', 'tpsfcube', 'bkgcube', 'tmodel', 'tonoff', 'tfit', 'tstat', 'tflux']

    hdr_sci = ''
    hdr_t = ''