import csv
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import os
from os import listdir
from os.path import isfile, isdir, join, expandvars

path = expandvars('$DATA')
folders = [f for f in listdir(path) if isdir(join(path, f))]

for i, folder in enumerate(folders):
    print(folder)
    tests = [f for f in listdir(join(path, folder)) if isfile(join(path, folder, f))]
    tests = sorted(tests)
    for j, test in enumerate(tests):
        data = pd.read_csv(join(path, folder, test), sep=' ')
        if j == 0:
            table = data
        else:
            table = table.append(data, sort=False)
    if i == 0:
        total = table
    else:
        total = total.append(table, sort=False)
table.to_csv(join(path, folder + '.txt'), index=False, header=True, sep=' ')
print(len(total))
total['offset'] = np.where(total['offset'] == 'gw', float(1.6), total['offset'])
#total.sort_index(axis=1, inplace=True)
total.to_csv(join(path, 'all.txt'), index=False, header=True, sep=' ')

print('exit')