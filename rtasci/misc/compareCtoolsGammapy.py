# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import pandas as pd
import seaborn as sns
from lib.RTAVisualise import *
from os.path import join, isdir, isfile
from os import listdir

datapath = join(os.path.expandvars('$DATA'), 'outputs/crab/timing')
try:
    isdir(datapath)
except ValueError:
    raise ValueError(f'Folder {datapath} not found.')

pngpath = join(datapath, 'png')
if not isdir(pngpath):
    print('Creating png folder...')
    os.mkdir(pngpath)

tables = [f for f in listdir(datapath) if isfile(join(datapath, f))]
print(f'Tables: {len(tables)}\n')

for i, table in enumerate(tables):
    print(f'Collecting data from {table}')
    data = pd.read_csv(join(datapath, table), sep=' ')
    data['tools'] = [table.replace('.csv', '') for n in range(len(data))]
    #print(f'5 rows {data[:5]}')
    if i == 0:
        print('start')
        total = data
    else:
        print('add')
        total = total.append(data, sort=False)
    print(f'data: {len(data)} and keys {len(data.keys())}')
    print(f'total: {len(total)} and keys {len(total.keys())}')


tools = total['tools'].drop_duplicates()
print(f'\n-----------\nanalysis: {len(tools)}')

for tool in tools:
    d1000 = total[total['texp'] == 1000]
    d1000 = d1000[d1000['tools'] == tool]
    d100 = total[total['texp'] == 100]
    d100 = d100[d100['tools'] == tool]
    d10 = total[total['texp'] == 10]
    d10 = d10[d10['tools'] == tool]


g = sns.FacetGrid(total, row="tools", col="ttotal", height=4, aspect=.5)
g.map(sns.barplot, "texp", "sqrt_ts")
