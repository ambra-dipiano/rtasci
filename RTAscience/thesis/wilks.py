# ----------------------------- #
# Copyright 2020 Ambra Di Piano #
# ----------------------------- # -------------------------------------------------- #
# Redistribution and use in source and binary forms, with or without modification,   #
# are permitted provided that the following conditions are met:                      #
# 1. Redistributions of source code must retain the above copyright notice,          #
# this list of conditions and the following disclaimer.                              #
# 2. Redistributions in binary form must reproduce the above copyright notice,       #
# this list of conditions and the following disclaimer in the documentation and/or   #
# other materials provided with the distribution.                                    #
# 3. Neither the name of the copyright holder nor the names of its contributors      #
# may be used to endorse or promote products derived from this software without      #
# specific prior written permission.                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. #
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,   #
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,     #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,      #
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE    #
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  #
# OF THE POSSIBILITY OF SUCH DAMAGE.                                                 #
# ---------------------------------------------------------------------------------- #

# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMpORTS ---!
from pkg_blindsearch import *
import csv
import os
import sys

# --------------------------------- SETUP --------------------------------- !!!

# initialize global count ---!
chunk = int(sys.argv[1])  # global count
trials = int(sys.argv[2])  # number of trials
count = int(sys.argv[3])  # starting count
# cpus ---!
nthreads = 1
os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
os.environ['MKL_NUM_THREADS'] = str(nthreads)

# ctools/cscripts parameters ---!
caldb = 'prod3b-v2'
irf = 'South_z40_0.5h'

texp = [1, 5, 10, 100]  # exposure times (s)
texp.sort()
tint = len(texp)
tmin = 0  # slewing time (s)
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 150.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 150.0  # selection maximum energy (TeV)
roi = 5  # region of interest (deg)

# conditions control ---!
checks = False  # prints checks info ---!
irf_degrade = False  # use degraded irf ---!
skip_exist = False  # if an output already exists it skips the step ---!
debug = False  # prints logfiles on terminal ---!
if_log = True  # saves logfiles ---!

# files ---!
model_pl = '/crab_z40.xml'
model_bkg = '/CTAIrfBackground.xml'
tcsv = 'time_slices.csv'
cfg = xmlConfig()
p = ConfigureXml(cfg)

# pointing ---!
pointing = (31.582, -53.211)  # pointing direction RA/DEC (deg)

# recap and dof ---!
dof, m2, m1 = getDof()
print('!!! *** !!! dof = ', m2, ' - ', m1, ' = ', dof)
print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! nominal prod:', caldb)
print('!!! *** !!! irf:', irf)
print('!!! *** !!! sim energy range: [', elow, ', ', ehigh, '] (TeV)')
print('!!! *** !!! selection energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi: ', roi, ' (deg)')
print('!!! *** !!! pointing:', pointing, ' (deg)')

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup model dof ---!i
if count == 0:
  model = p.getWorkingDir() + model_pl
  mObj = ManageXml(model)
  mObj.prmsFreeFix()
  mObj.closeXml()
# setup trials obj ---!
tObj = Analysis()
tObj.nthreads = nthreads
tObj.pointing = pointing
tObj.roi = roi
tObj.e = [elow, ehigh]
tObj.tmax = max(tmax)
tObj.debug = debug
tObj.if_log = if_log
# degrade IRF if required ---!
if irf_degrade:
  tObj.caldb = caldb.replace('prod', 'degr')
else:
  tObj.caldb = caldb
tObj.irf = irf
print('!!! check ---- caldb:', tObj.caldb) if checks is True else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
  count += 1
  tObj.seed = count
  print('\n\n!!! ************ STARTING TRIAL %d ************ !!!' %count) if checks is True else None
  print('!!! check ---- seed=', tObj.seed) if checks is True else None

  # --------------------------------- CHECK SKIP --------------------------------- !!!

  ID = 'ID%06d' % count
  csvName = p.getCsvDir() + 'bkg_%ds_chunk%02d.csv' % (texp[-1], chunk)
  if os.path.isfile(csvName) and skip_exist:
    skip = checkTrialId(csvName, ID)
  else:
    skip = False
  if skip_exist is True and skip is True:
    continue

  # --------------------------------- SIMULATION --------------------------------- !!!

  # attach ID to fileroot ---!
  f = 'bkg%06d' % (count)
  if irf_degrade:
    f += 'irf'
  # simulate ---!
  model = p.getWorkingDir() + model_bkg
  tObj.model = model
  event = p.getSimDir() + f + ".fits"
  if os.path.isfile(event):
    os.remove(event)
  tObj.output = event
  tObj.eventSim()
  print('!!! check ---- simulation=', event) if checks is True else None

  # --------------------------------- 2° LOOP :: texp --------------------------------- !!!

  # --------------------------------- SELECTION --------------------------------- !!!

  tObj.e = [emin, emax]
  for i in range(tint):
    print('\n\n!!! ************ STARTING TEXP %d ************ !!!\n\n' % texp[i]) if checks is True else None
    tObj.t = [tmin, tmax[i]]
    event_selected = event.replace(p.getSimDir(), p.getSelectDir()).replace('bkg%06d' %count,
                                                                            'texp%ds_bkg%06d' %(texp[i], count))
    prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
    if os.path.isfile(event_selected):
      os.remove(event_selected)
    tObj.input = event
    tObj.output = event_selected
    tObj.eventSelect(prefix=prefix)
    print('!!! check ---- selection: ', event_selected) if checks is True else None

    # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

    model = p.getWorkingDir() + model_pl
    likeXml = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.fits', '_like.xml')
    if os.path.isfile(likeXml):
      os.remove(likeXml)
    tObj.input = event_selected
    tObj.model = model
    tObj.output = likeXml
    tObj.maxLikelihood()
    likeObj = ManageXml(likeXml)
    print('!!! check ---- max likelihood: ', likeXml) if checks is True else None

    # --------------------------------- BEST FIT TSV --------------------------------- !!!

    ts_list, ts = ([] for j in range(2))
    ts_list.append(likeObj.loadTs())

    # only first elem ---!
    ts.append(ts_list[0][0])

    # --------------------------------- BEST FIT RA & DEC --------------------------------- !!!

    ra_list, ra_fit, dec_list, dec_fit = ([] for j in range(4))
    coord = likeObj.loadRaDec()
    ra_list.append(coord[0])
    dec_list.append(coord[1])

    # only first elem ---!
    ra_fit.append(ra_list[0][0])
    dec_fit.append(dec_list[0][0])

    # --------------------------------- BEST FIT SPECTRAL --------------------------------- !!!

    pref_list, pref, index_list, index, pivot_list, pivot = ([] for j in range(6))
    spectral = likeObj.loadSpectral()
    index_list.append(spectral[0])
    pref_list.append(spectral[1])
    pivot_list.append(spectral[2])

    # only first elem ---!
    index.append(index_list[0][0])
    pref.append(pref_list[0][0])
    pivot.append(pivot_list[0][0])

    # --------------------------------- INTEGRATED FLUX --------------------------------- !!!

    flux_ph = []
    flux_ph.append(tObj.photonFluxPowerLaw(index[0], pref[0], pivot[0]))  # E (MeV)

    # --------------------------------- CLOSE LIKE XML --------------------------------- !!!

    likeObj.closeXml()

    # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

    header = '#trial,RA_fit,DEC_fit,flux_ph,TS\n'
    ID = 'ID%06d' % count
    csvName = p.getCsvDir() + 'bkg_%ds_chunk%02d.csv' % (texp[i], chunk)

    row = []
    print('\n\n!!! ---------- check trial:', count) if checks is True else None
    print('!!! ----- check texp:', texp[i]) if checks is True else None
    print('!!! *** check ts:', ts[0][0]) if checks is True else None

    row.append([ID, ra_fit[0], dec_fit[0], flux_ph[0], ts[0]])
    print('!!! check row: seed %d --- texp' %i, texp[i], 's =====', row) if checks is True else None
    if os.path.isfile(csvName):
      with open(csvName, 'a') as f:
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    else:
      with open(csvName, 'w+') as f:
        f.write(header)
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    print('!!! check ---- data file: ', csvName) if checks is True else None

  # --------------------------------- CLEAR SPACE --------------------------------- !!!

  print('!!! check ---- ', count, ') trial done...') if checks is True else None
  if int(count) != 1:
    os.system('rm ' + p.getSimDir() + '*bkg%06d*' % count)
    os.system('rm ' + p.getSelectDir() + '*bkg%06d*' % count)
    os.system('rm ' + p.getDetDir() + '*bkg%06d*' % count)

print('\n\n!!! ================== END ================== !!!\n\n')


