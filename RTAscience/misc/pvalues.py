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
from os.path import join, expandvars
from RTAscience.lib.RTAStats import *

path = expandvars('$DATA/outputs/LEO')
png_path = join(path, 'png/')
if not os.path.isdir(png_path):
  os.mkdir(png_path)


filename = 'chi2sample_s1e7_r0-36_n100.csv'
# load DataFrame and column names ---!
df = pd.read_csv(join(path, filename))
cols = list(df.columns)
trials = len(df[cols[0]])
# drop duplicates ---!
df.sort_values(cols[0], inplace=True)
# set arrays ---!
values = np.array(df[cols[0]])
values.sort()

nbin = 100
print(f'min = {min(values)}\nmax = {max(values)}\nlen = {len(values)}')

print(len(values[values > 36]), len(values[values >= 36]))

print(values[values >= 36])

# -------------------------------- PLOT ---!

fig, ax = ts_wilks(values, len(values), nbin=nbin, figsize=(7, 8), xrange=(0,36), title='distribution', show=False, usetex=False, filename=png_path + filename.replace('.csv', '_wilks.png'), overlay=None, write_data=True)

fig, ax = p_values(values, len(values), nbin=nbin, figsize=(7, 8), xrange=(0,36), title='pvalues', show=False, usetex=False, filename=png_path + filename.replace('.csv', '_pvalues.png'), overlay=None, sigma5=False, write_data=True)

#fig, ax = ts_wilks_cumulative(values, len(values), nbin=nbin, figsize=(7, 8), xrange=(0,0.3), show=False, usetex=False, title='cumulative', filename=png_path + filename.replace('.csv', '_cumulative.png'), overlay=None)
