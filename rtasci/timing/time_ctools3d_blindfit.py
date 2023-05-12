# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import time
import sys
texp = sys.argv[1]
first = sys.argv[2]

# start timing
t = time.time()
clock0 = time.time()
import numpy as np
from rtasci.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from rtasci.lib.RTAUtils import phflux_powerlaw
from rtasci.lib.RTAManageXml import ManageXml
timport = time.time() - t
print(f'Imports : {timport} s\n')

t = time.time()
rootpath = "/data01/homes/cta/gammapy_integration"
obspath = f'{rootpath}/DATA/obs/crab/'
rtapath = f'{rootpath}/DATA/rta_products/crab/'
filename = f'{obspath}crab_offax_texp{texp}s_n01.fits'
skyname = filename.replace(obspath,rtapath).replace('.fits', '_skymap.fits')
detname = skyname.replace('_skymap.fits',f'_model.xml')
fitname = detname.replace('_model.xml','_fit.xml')
print(f'Fits: {filename.replace(obspath, "")}\n')
tsetup = time.time() - t
print(f'Setup : {tsetup} s\n')

# initialise + skymap
t = time.time()
analysis = RTACtoolsAnalysis()
analysis.nthreads = 1
analysis.caldb = 'prod3b-v2'
analysis.irf = 'South_z20_0.5h'
analysis.e = [0.05, 20]
analysis.t = [0, texp]
analysis.roi = 5
analysis.input = filename
analysis.output = skyname
analysis.run_skymap()
tsky = time.time() - t
print(f'Skymap: {tsky} s\n')

# initialise
t = time.time()
analysis.sigma = 3
analysis.max_src = 1
analysis.input = skyname
analysis.output = detname
analysis.run_blindsearch()
tblind = time.time() - t
print(f'Blind-search: {tblind} s\n')

# get candidate and modify model
t = time.time()
detection = ManageXml(detname)
detection.modXml(overwrite=True)
detection.setTsTrue() 
detection.parametersFreeFixed(src_free=['Prefactor'])
detection.closeXml()
tmodel = time.time() - t
print(f'Modelling: {tmodel} s\n')

# fitting
t = time.time()
analysis.input = filename
analysis.model = detname
analysis.output = fitname
analysis.run_maxlikelihood()
tfit = time.time() - t
print(f'Fitting: {tfit} s\n')

# statistics
t = time.time()
results = ManageXml(fitname)
try:
    coords = results.getRaDec()
    ra = coords[0][0]
    dec = coords[1][0]
    ts = results.getTs()[0]
except IndexError:
    raise Warning('No candidates found.')
print(f'Hotspots:{coords}\n')
print(f'sqrt_ts: {np.sqrt(ts)}')
tstat = time.time() - t
print(f'Statistics: {tstat} s\n')

# flux
t = time.time()
spectra = results.getSpectral()
index, pref, pivot = spectra[0][0], spectra[1][0], spectra[2][0]
err = results.getPrefError()[0]
phflux = phflux_powerlaw(index, pref, pivot, analysis.e, unit='TeV')
phflux_err = phflux_powerlaw(index, err, pivot, analysis.e, unit='TeV')
print(f'PH-FLUX {phflux} +/- {phflux_err}\n')
tflux = time.time() - t
print(f'Flux points : {tflux} s\n')

ttotal = time.time() - clock0
print(f'Total time: {ttotal} s\n')
print('\n\n-----------------------------------------------------\n\n')

logname = f'{rootpath}/DATA/outputs/crab/ctools3d_blindfit.csv'
row = f'{texp} {np.sqrt(ts)} {phflux} {phflux_err} {ra} {dec} {ttotal} {timport} {tsetup} {tsky} {tblind} {tmodel} {tfit} {tstat} {tflux}\n'
if first == 'True':
    hdr = 'texp sqrt_ts flux flux_err ra dec ttotal timport tsetup tsky tblind tmodel tfit tstat tflux\n'
    log = open(logname, 'w+')
    log.write(hdr)
    log.write(row)
    log.close()
else:
    log = open(logname, 'a')
    log.write(row)
    log.close()

print(row)