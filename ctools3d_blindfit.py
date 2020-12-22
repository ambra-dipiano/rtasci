# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import time
t = time.time()
clock0 = time.time()
import numpy as np
from lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from lib.RTAUtils import phflux_powerlaw
from lib.RTAManageXml import ManageXml
print(f'Imports : {time.time() - t} s\n')

t = time.time()
obspath = '/home/ambra/Desktop/CTA/projects/DATA/obs/crab'
rtapath = '/home/ambra/Desktop/CTA/projects/DATA/rta_products/crab'
filename = f'{obspath}/crab_onax.fits'
print(f'Fits: {filename.replace(obspath, "")}\n')
skyname = filename.replace(obspath,rtapath).replace('.fits', '_skymap.fits')
detname = skyname.replace('_skymap.fits',f'_model.xml')
fitname = detname.replace('_model.xml','_fit.xml')
print(f'Setup : {time.time() - t} s\n')

# initialise + skymap
t = time.time()
analysis = RTACtoolsAnalysis()
analysis.nthreads = 1
analysis.caldb = 'prod3b-v2'
analysis.irf = 'South_z20_0.5h'
analysis.e = [0.05, 20]
analysis.input = filename
analysis.output = skyname
analysis.run_skymap()
print(f'Skymap: {time.time() - t} s\n')

# initialise
t = time.time()
analysis.sigma = 3
analysis.max_src = 1
analysis.input = skyname
analysis.output = detname
analysis.run_blindsearch()
print(f'Blind-search: {time.time() - t} s\n')

# get candidate and modify model
t = time.time()
detection = ManageXml(detname)
detection.modXml(overwrite=True)
detection.setTsTrue() 
detection.parametersFreeFixed(src_free=['Prefactor'])
detection.closeXml()
print(f'Modelling: {time.time() - t} s\n')

# fitting
t = time.time()
analysis.input = filename
analysis.model = detname
analysis.output = fitname
analysis.run_maxlikelihood()
print(f'Fitting: {time.time() - t} s\n')

# statistics
t = time.time()
results = ManageXml(fitname)
try:
    coords = results.getRaDec()
    ts = results.getTs()[0]
except IndexError:
    raise Warning('No candidates found.')
print(f'Hotspots:{coords}\n')
print(f'sqrt_ts: {np.sqrt(ts)}')
print(f'Statistics: {time.time() - t} s\n')

# flux
t = time.time()
spectra = results.getSpectral()
index, pref, pivot = spectra[0][0], spectra[1][0], spectra[2][0]
err = results.getPrefError()[0]
phflux = phflux_powerlaw(index, pref, pivot, analysis.e, unit='TeV')
phflux_err = phflux_powerlaw(index, err, pivot, analysis.e, unit='TeV')
print(f'PH-FLUX {phflux} +/- {phflux_err}\n')
print(f'Flux points : {time.time() - t} s\n')
#print(result)

print(f'Total time: {time.time() - clock0} s\n')

