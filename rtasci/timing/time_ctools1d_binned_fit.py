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
from rtasci.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from rtasci.lib.RTAUtils import *
from rtasci.lib.RTAManageXml import ManageXml
timport = time.time() - t
print(f'Imports : {timport} s\n')

t = time.time()
edisp = True
rootpath = "/data01/homes/cta/gammapy_integration"
obspath = f'{rootpath}/DATA/obs/crab/'
rtapath = f'{rootpath}/DATA/rta_products/crab/'
modelpath = f'{rootpath}/DATA/models/'
filename = f'{obspath}crab_offax_texp{texp}s_n01.fits'
model = f'{modelpath}crab.xml'
onoff_obs = filename.replace(obspath,rtapath).replace('.fits','_cspha.xml')
onoff_model = onoff_obs.replace('.xml','_model.xml')
fitname = onoff_obs.replace('.xml','_fit.xml')
print(f'Fits: {filename.replace(obspath, "")}\n')
tsetup = time.time() - t
print(f'Setup : {tsetup} s\n')

# set model
t = time.time()
xml = ManageXml(model)
xml.setTsTrue() 
xml.parametersFreeFixed(src_free=['Prefactor'])
target = xml.getRaDec()
xml.closeXml()
tmodel = time.time() - t
print(f'Modelling: {tmodel} s\n')

# onoff
t = time.time()
analysis = RTACtoolsAnalysis()
analysis.nthreads = 1
analysis.caldb = 'prod3b-v2'
analysis.irf = 'South_z20_0.5h'
analysis.e = [0.05, 20]
analysis.t = [0, texp]
analysis.roi = 5
analysis.input = filename
analysis.model = model
analysis.src_name = 'Crab'
analysis.target = [target[0][0], target[1][0]]
analysis.output = onoff_obs
analysis.run_onoff(prefix=onoff_obs.replace('.xml',''), ebins=30, etruemin=0.03, etruemax=30, etruebins=40, ebins_alg='LOG', maxoffset=2.5)
tonoff = time.time() - t
print(f'Onoff: {tonoff} s\n')

# onoff
t = time.time()
analysis.input = onoff_obs
analysis.model = onoff_model
analysis.output = fitname
analysis.run_maxlikelihood(edisp=edisp)
tfit = time.time() - t
print(f'Fitting: {tfit} s\n')

# statistics
t = time.time()
results = ManageXml(fitname)
try:
    ts = results.getTs()[0]
except IndexError:
    raise Warning('Target not found.')
print(f'sqrt_ts: {np.sqrt(ts)}')
tstat = time.time() - t
print(f'Statistics: {tstat} s\n')

# flux
t = time.time()
spectra = results.getSpectral()
index, pref, pivot = spectra[0][0], spectra[1][0], spectra[2][0]
err = results.getPrefError()[0]
print(analysis.e)
phflux = phflux_powerlaw(index, pref, pivot, analysis.e, unit='TeV')
phflux_err = phflux_powerlaw(index, err, pivot, analysis.e, unit='TeV')
print(f'PH-FLUX {phflux} +/- {phflux_err}\n')
tflux = time.time() - t
print(f'Flux points : {tflux} s\n')

ttotal = time.time() - clock0
print(f'Total time: {ttotal} s\n')
print('\n\n-----------------------------------------------------\n\n')

logname = f'{rootpath}/DATA/outputs/crab/ctools1d_binned_fit.csv'
row = f'{texp} {np.sqrt(ts)} {phflux} {phflux_err} {ttotal} {timport} {tsetup} {tmodel} {tonoff} {tfit} {tstat} {tflux}\n'
if first == 'True':
    hdr = 'texp sqrt_ts flux flux_err ttotal timport tsetup tmodel tonoff tfit tstat tflux\n'
    log = open(logname, 'w+')
    log.write(hdr)
    log.write(row)
    log.close()
else:
    log = open(logname, 'a')
    log.write(row)
    log.close()

print(row)