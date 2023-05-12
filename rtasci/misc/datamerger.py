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
from os import listdir
from os.path import isfile, isdir, join, expandvars

path = expandvars('$DATA')
folders = [f for f in listdir(path) if isdir(join(path, f))]

for i, folder in enumerate(folders):
    print(folder)
    tests = [f for f in listdir(join(path, folder)) if isfile(join(path, folder, f)) and '.txt' in f]
    tests = sorted(tests)
    for j, test in enumerate(tests):
        data = pd.read_csv(join(path, folder, test), sep=' ')
        if j == 0:
            table = data
        else:
            table = table.append(data, sort=False)
    print(len(table))
    table.to_csv(join(path, folder + '.txt'), index=False, header=True, sep=' ')
    if i == 0:
        total = table
    else:
        total = total.append(table, sort=False)
print(len(total))
total['offset'] = np.where(total['offset'] == 'gw', float(1.6), total['offset'])
#total.sort_index(axis=1, inplace=True)
total.to_csv(join(path, 'all.txt'), index=False, header=True, sep=' ', na_rep=np.nan)

print('exit')