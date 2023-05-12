# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from os import listdir
from os.path import isfile, isdir, join, expandvars

# target
trueRA = 33.057
trueDEC = -51.841
true_coord = SkyCoord(ra = trueRA*u.deg, dec = trueDEC*u.deg, frame='fk5')

# data
path = expandvars('$DATA')
folders = [f for f in listdir(path) if isdir(join(path, f)) and 'offGW' in f and 'del50' in f]
folders = sorted(folders)
for i, folder in enumerate(folders):
    print(folder)
    tests = [f for f in listdir(join(path, folder)) if isfile(join(path, folder, f)) and '.csv' in f]
    tests = sorted(tests)
    for j, test in enumerate(tests):
        data = pd.read_csv(join(path, folder, test), sep=',')
        data['scaleflux'] = i+1
        data = data.rename(columns={'#trial': 'runid', 'RA_det': 'ra', 'DEC_det': 'dec', 'flux_ph': 'flux', 'TS': 'sqrt_ts'})
        data['seed'] = [int(seed.replace('ID','')) for seed in data['runid']]
        data['sqrt_ts'] = [np.sqrt(float(ts)) for ts in np.array(data['sqrt_ts'])]
        data['runid'] = 'run0406_ID000126'
        data['offset'] = 1.638
        data['delay'] = 50
        data['caldb'] = 'degr3b-v2'
        data['irf'] = 'South_z40_0.5h'
        data = data.drop(['Ndet', 'Nsrc', 'RA_fit', 'DEC_fit', 'sigma'], axis=1)
        dist = []
        for k in range(len(data)):
            dist.append(float(true_coord.separation(SkyCoord(ra=np.array(data['ra'])[k], dec=np.array(data['dec'])[k], unit='deg', frame='fk5')).deg))
        print(len(dist), max(dist), min(dist))
        data['dist'] = dist
        #data.sort_index(axis=1, inplace=True)
        columnsTitles = ['runid', 'seed', 'texp', 'sqrt_ts', 'flux', 'ra', 'dec', 'dist', 'offset', 'delay', 'scaleflux', 'caldb', 'irf']
        data = data.reindex(columns=columnsTitles)
        data.to_csv(join(path, folder, test.replace('.csv', '.txt')), sep=' ', index=False, header=True, na_rep='nan')
        if j == 0:
            table = data
        else:
            table = table.append(data, sort=False)
    table.to_csv(join(path, folder + '.txt'), sep=' ', index=False, header=True, na_rep='nan')
    del table