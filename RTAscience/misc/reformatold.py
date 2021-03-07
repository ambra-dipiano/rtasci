import csv
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import os
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
    tests = [f for f in listdir(join(path, folder)) if isfile(join(path, folder, f))]
    tests = sorted(tests)
    for j, test in enumerate(tests):
        data = pd.read_csv(join(path, folder, test), sep=',')
        data['scaleflux'] = i+1
        if i == 0:
            if j == 0:
                table = data
            else:
                table = table.append(data, sort=False)
        else:
            table = table.append(data, sort=False)
    table = table.rename(columns={'#trial': 'runid', 'TS': 'sqrt_ts', 'RA_det': 'ra', 'DEC_det': 'dec', 'flux_ph': 'flux'})
    table['seed'] = [int(seed.replace('ID','')) for seed in table['runid']]
    table['sqrt_ts'] = [np.sqrt(ts) for ts in table['sqrt_ts']]
    table['runid'] = 'run0406_ID000126'
    table['offset'] = 2
    table['delay'] = 50
    table['caldb'] = 'degr3b-v2'
    table['irf'] = 'South_z40_0.5h'
    table = table.drop(['Ndet', 'Nsrc', 'RA_fit', 'DEC_fit', 'sigma'], axis=1)
    dist = []
    for k in range(len(table)):
        print(k)
        dist.append(float(true_coord.separation(SkyCoord(ra=np.array(table['ra'])[k], dec=np.array(table['dec'])[k], unit='deg', frame='fk5')).deg))
    print(len(dist), max(dist), min(dist))
    table['dist'] = dist
    table = table.sort('seed')
    #table.sort_index(axis=1, inplace=True)
    columnsTitles = ['runid', 'seed', 'texp', 'sqrt_ts', 'flux', 'ra', 'dec', 'dist', 'offset', 'delay', 'scaleflux', 'caldb', 'irf']
    table = table.reindex(columns=columnsTitles)
    table.to_csv(join(path, folder + '.txt'), sep=' ', index=False, header=True)