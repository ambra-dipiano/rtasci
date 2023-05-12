# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
from rtasci.lib.RTAStats import *

dof = 1
folder = 'tesi_bkg_1e6_nominal_%ddof_fullE_old' % dof
root = '/home/ambra/Desktop/CTA/projects/rta-pipe'
path = root + '/archive_tests/tesi_02_bkg/' + folder + '/run0406_bkg/run0406_ID000126/csv/'
png_path = root + '/archive_tests/paper2021/png/'
if not os.path.isdir(png_path):
  os.mkdir(png_path)

Nchunk = 20

texp = [1, 5, 10, 100]
sigma = [5]
chunk = [i + 1 for i in range(Nchunk)]

# csvName[texp][chunk]
csvName = [[] * i for i in range(len(texp))]
for i in range(len(chunk)):
  for j in range(len(texp)):
    csvName[j].append('bkg_%ds_chunk%02d.csv' % (texp[j], chunk[i]))

# merge files ---!
csvMerged = []
for j in range(len(texp)):
  csvMerged.append('bkg_%ds.csv' % texp[j])

  fout = open(path + csvMerged[j], 'w+')
  # first file ---!
  for line in open(path + csvName[j][0]):
    fout.write(line)
  # remaining files ---!
  for i in range(len(chunk) - 1):
    f = open(path + csvName[j][i + 1])
    next(f)  # skip the header ---!
    for line in f:
      fout.write(line)
    f.close()
  fout.close()

print('data files merge completed')

show = False
fontsize = 18
for n in range(len(texp)):
  filename = csvMerged[n]
  print('!======== texp = ', texp[n], '(s) ========!')
  # load DataFrame and column names ---!
  df = pd.read_csv(path + filename)
  cols = list(df.columns)
  trials = len(df[cols[0]])
  print('verify trials = ', trials)
  # drop duplicates ---!
  df.sort_values(cols[0], inplace=True)
  # dropping ALL duplicte values
  df.drop_duplicates(subset=cols[0], keep='last', inplace=True)
  trials = len(df[cols[0]])
  print('verify trials = ', trials)
  # drop NaN ---!
  # df = df.dropna()

  # set arrays ---!
  trial = np.array(df[cols[0]])
  tsv = np.array(df[cols[-1]])

  tsv.sort()

  wbin = 1
  nbin = int(tsv.max() / wbin)
  if nbin == 0:
    nbin = 1
  print('ts bin:', nbin)

  # -------------------------------- STATS ---!

  ts = []
  for i in range(trials):
    if tsv[i] < 0.0 or tsv[i] == np.nan:
      ts.append(0.0)
    else:
      ts.append(tsv[i])

  # chi2, chi2r = chi2_reduced(ts, trials, df=dof, nbin=nbin, width=wbin, var=False)
  # print('var=False; chi2=', chi2, '; chi2r=', chi2r)


  # -------------------------------- PLOT ---!

  fig, ax = ts_wilks(ts, len(ts), df=dof, nbin=nbin, width=wbin, figsize=(10, 6), fontsize=fontsize,
                     title='prod3b-v2: South_z40_0.5h (texp=%ds)' % texp[n], show=True, usetex=False,
                     filename=png_path + filename.replace('.csv', '_wilks.png'), ylim=(1e-7, 2e0), 
                     xlim=(0.0, 30))

  fig, ax = p_values(ts, len(ts), df=dof, nbin=nbin, width=wbin, figsize=(10, 6), fontsize=fontsize,
                     title='prod3b-v2: South_z40_0.5h (texp=%ds)' % texp[n], show=False, usetex=False,
                     filename=png_path + filename.replace('.csv', '_pvalues.png'), ylim=(1e-7, 2e0), 
                     xlim=(0.0, 30))

  fig, ax = ts_wilks_cumulative(ts, len(ts), df=dof, nbin=nbin, width=wbin, figsize=(10, 6),
                                fontsize=fontsize, show=False, usetex=False,
                                title='prod3b-v2: South_z40_0.5h (texp=%ds)' % texp[n],
                                filename=png_path + filename.replace('.csv', '_cumulative.png'))



