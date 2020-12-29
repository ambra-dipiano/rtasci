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
from lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from lib.RTAUtils import *
from lib.RTAManageXml import ManageXml
print(f'Imports : {time.time() - t} s\n')

t = time.time()
obspath = '/home/ambra/Desktop/CTA/projects/DATA/obs/crab'
rtapath = '/home/ambra/Desktop/CTA/projects/DATA/rta_products/crab'
modelpath = '/home/ambra/Desktop/CTA/projects/DATA/models'
filename = f'{obspath}/crab_offax.fits'
model = f'{modelpath}/crab.xml'
onoff_obs = filename.replace(obspath,rtapath).replace('.fits','_cspha.xml')
onoff_model = onoff_obs.replace('.xml','_model.xml')
fitname = onoff_obs.replace('.xml','_fit.xml')
print(f'Fits: {filename.replace(obspath, "")}\n')
print(f'Setup : {time.time() - t} s\n')

# set model
t = time.time()
xml = ManageXml(model)
xml.setTsTrue() 
xml.parametersFreeFixed(src_free=['Prefactor'])
target = xml.getRaDec()
xml.closeXml()
print(f'Modelling: {time.time() - t} s\n')

# onoff
t = time.time()
analysis = RTACtoolsAnalysis()
analysis.nthreads = 1
analysis.caldb = 'prod3b-v2'
analysis.irf = 'South_z20_0.5h'
analysis.e = [0.05, 20]
analysis.input = filename
analysis.model = model
analysis.src_name = 'Crab'
analysis.target = [target[0][0], target[1][0]]
analysis.output = onoff_obs
analysis.run_onoff()
print(f'Onoff: {time.time() - t} s\n')

# onoff
t = time.time()
analysis.input = onoff_obs
analysis.model = onoff_model
analysis.output = fitname
analysis.run_maxlikelihood()
print(f'Fitting: {time.time() - t} s\n')

# statistics
t = time.time()
results = ManageXml(fitname)
try:
    ts = results.getTs()[0]
except IndexError:
    raise Warning('Target not found.')
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


print(f'Total time: {time.time() - clock0} s\n')
