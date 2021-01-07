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
from lib.RTACtoolsAnalysis import RTACtoolsAnalysis, make_obslist
from lib.RTAUtils import *
from lib.RTAManageXml import ManageXml
timport = time.time() - t
print(f'Imports : {timport} s\n')

t = time.time()
obspath = '/home/ambra/Desktop/CTA/projects/DATA/selections/crab/'
rtapath = '/home/ambra/Desktop/CTA/projects/DATA/rta_products/crab/'
modelpath = '/home/ambra/Desktop/CTA/projects/DATA/models/'
filename = f'{obspath}crab_offax_texp{texp}s_n01.fits'
obslist = filename.replace('.fits', '_obs.xml')
cube = filename.replace(obspath,rtapath).replace('.fits','_cube.xml')
expcube = cube.replace('.xml', '_exp.xml')
psfcube = cube.replace('.xml', '_psf.xml')
bkgcube = cube.replace('.xml', '_bkg.xml')
inmodel = f'{modelpath}crab.xml'
outmodel = cube.replace('.xml', '_model.xml')
fitname = filename.replace(obspath,rtapath).replace('.fits', '_fit.xml')
print(f'Fits: {filename.replace(obspath, "")}\n')
tsetup = time.time() - t
print(f'Setup : {tsetup} s\n')

t = time.time()
make_obslist(obslist, filename, 'Crab')
tobs = time.time() - t
print(f'Observation: {tobs} s\n')

# initialise + binning
t = time.time()
analysis = RTACtoolsAnalysis()
analysis.nthreads = 1
analysis.caldb = 'prod3b-v2'
analysis.irf = 'South_z20_0.5h'
analysis.e = [0.05, 20]
analysis.usepnt = True
analysis.input = obslist
analysis.output = cube
analysis.run_binning(prefix=cube.replace('.xml','_'), ebins=10)
tcube = time.time() - t
print(f'Binning: {tcube} s\n')

# exp cube
t = time.time()
analysis.input = obslist
analysis.output = expcube
analysis.run_expcube(cube=cube.replace('.xml','_cta.fits'), ebins=10)
texpcube = time.time() - t
print(f'Exp. cube: {texpcube} s\n')

# psf cube
t = time.time()
analysis.input = obslist
analysis.output = psfcube
analysis.run_psfcube(cube=cube.replace('.xml','_cta.fits'), ebins=10)
tpsfcube = time.time() - t
print(f'Psf. cube: {tpsfcube} s\n')

# bkg cube
t = time.time()
analysis.input = obslist
analysis.model = inmodel
analysis.output = bkgcube
analysis.run_bkgcube(cube=cube.replace('.xml','_cta.fits'), model=outmodel)
tbkgcube = time.time() - t
print(f'Bkg. cube: {tbkgcube} s\n')

# set model
t = time.time()
xml = ManageXml(outmodel)
xml.setTsTrue() 
xml.parametersFreeFixed(src_free=['Prefactor'])
xml.closeXml()
tmodel = time.time() - t
print(f'Modelling: {tmodel} s\n')

# fitting
analysis.input = cube
analysis.model = outmodel
analysis.output = fitname
analysis.run_maxlikelihood(binned=True, exp=expcube, psf=psfcube, bkg=bkgcube)
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
phflux = phflux_powerlaw(index, pref, pivot, analysis.e, unit='TeV')
phflux_err = phflux_powerlaw(index, err, pivot, analysis.e, unit='TeV')
print(f'PH-FLUX {phflux} +/- {phflux_err}\n')
tflux = time.time() - t
print(f'Flux points : {tflux} s\n')

ttotal = time.time() - clock0
print(f'Total time: {ttotal} s\n')
print('\n\n-----------------------------------------------------\n\n')

logname = f'/home/ambra/Desktop/CTA/projects/DATA/outputs/crab/ctools3d_binned_fit.csv'
if first:
    hdr = 'texp sqrt_ts flux flux_err ttotal timport tsetup tobs tcube texpcube tpsfcube tbkgcube tmodel tfit tstat tflux\n'
    log = open(logname, 'w+')
    log.write(hdr)
    log.write(f'{texp} {np.sqrt(ts)} {phflux} {phlux_err} {ttotal} {timport} {tsetup} {tmodel} {tonoff} {tfit} {tstat} {tflux}\n')
    log.close()
else:
    log = open(logname, 'a')
    log.write(f'{texp} {np.sqrt(ts)} {phflux} {phflux_err} {ttotal} {timport} {tsetup} {tobs} {tcube} {texpcube} {tepsfcube} {tbkgcube} {tmodel} {tfit} {tstat} {tflux}\n')
    log.close()



