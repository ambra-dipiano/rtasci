# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
from os.path import join, expandvars
from RTAscience.lib.RTAStats import *

path = expandvars('$DATA/outputs/LEO')
png_path = join(path, 'png/')
if not os.path.isdir(png_path):
  os.mkdir(png_path)


filename = 'reconstruction_errors_2.csv'
os.system(f'head {filename}')
# load DataFrame and column names ---!
df = pd.read_csv(join(path, filename))
cols = list(df.columns)
trials = len(df[cols[0]])
print('verify trials = ', trials)
# drop duplicates ---!
df.sort_values(cols[0], inplace=True)
# set arrays ---!
values = np.array(df[cols[0]])
values.sort()

values = np.genfromtxt(join(path, filename), usecols=(0), skip_header=0, dtype=float)

nbin = 100

# -------------------------------- PLOT ---!

fig, ax = ts_wilks(values, len(values), nbin=nbin, figsize=(10, 6), title='distribution', show=False, usetex=False, filename=png_path + filename.replace('.csv', '_wilks.png'))

fig, ax = p_values(values, len(values), nbin=nbin, figsize=(10, 6), title='pvalues', show=False, usetex=False, filename=png_path + filename.replace('.csv', '_pvalues.png'))

fig, ax = ts_wilks_cumulative(values, len(values), nbin=nbin, figsize=(10, 6), show=False, usetex=False, title='cumulative', filename=png_path + filename.replace('.csv', '_cumulative.png'))
