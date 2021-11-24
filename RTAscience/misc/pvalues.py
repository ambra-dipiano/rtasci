# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
from os.path import join
from RTAscience.lib.RTAStats import *

path = '/data01/homes/cta/gammapy_integration/DATA/outputs/LEO'
png_path = join(path, 'png/')
if not os.path.isdir(png_path):
  os.mkdir(png_path)


filename = 'reconstruction_errors.csv'
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

nbin = 40
wbin = values.max() / nbin
if nbin == 0:
    nbin = 1
print('width bin:', wbin, 'number of bins:', nbin)

# -------------------------------- PLOT ---!

fig, ax = ts_wilks(values, len(values), nbin=nbin, width=wbin, figsize=(10, 6), title='distribution', show=False, usetex=False, filename=png_path + filename.replace('.csv', '_wilks.png'))

fig, ax = p_values(values, len(values), nbin=nbin, width=wbin, figsize=(10, 6), title='pvalues', show=False, usetex=False, filename=png_path + filename.replace('.csv', '_pvalues.png'))

#fig, ax = ts_wilks_cumulative(values, len(values), nbin=nbin, width=wbin, figsize=(10, 6), show=False, usetex=False, title='cumulative', filename=png_path + filename.replace('.csv', '_cumulative.png'))
