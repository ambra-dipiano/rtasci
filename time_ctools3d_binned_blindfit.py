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
import os
texp = sys.argv[1]
first = sys.argv[2]

# start timing
t = time.time()
clock0 = time.time()
import numpy as np
from lib.RTACtoolsAnalysis import RTACtoolsAnalysis, make_obslist
from lib.RTAUtils import phflux_powerlaw
from lib.RTAManageXml import ManageXml
timport = time.time() - t
print(f'Imports : {timport} s\n')

t = time.time()
rootpath = str(os.path.dirname(os.path.abspath(__file__))).replace('cta-sag-sci', '')
obspath = f'{rootpath}/DATA/selections/crab/'
rtapath = f'{rootpath}/DATA/rta_products/crab/'
filename = f'{obspath}crab_offax_texp{texp}s_n01.fits'
skyname = filename.replace(obspath,rtapath).replace('.fits', '_skymap.fits')
detname = skyname.replace('_skymap.fits',f'_model.xml')
obslist = filename.replace('.fits', '_obs.xml')
cube = filename.replace(obspath,rtapath).replace('.fits','_cube.xml')
expcube = cube.replace('.xml', '_exp.xml')
psfcube = cube.replace('.xml', '_psf.xml')
bkgcube = cube.replace('.xml', '_bkg.xml')
outmodel = cube.replace('.xml', '_model.xml')
fitname = filename.replace(obspath,rtapath).replace('.fits', '_fit.xml')
print(f'Fits: {filename.replace(obspath, "")}\n')
tsetup = time.time() - t
print(f'Setup : {tsetup} s\n')

# observaition
t = time.time()
make_obslist(obslist, items=filename, names='Crab', instruments='CTA')
tobs = time.time() - t
print(f'Observation: {tobs} s\n')

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
detection.setInstrument()
detection.parametersFreeFixed(src_free=['Prefactor'])
try:
    coords = detection.getRaDec()
    ra = coords[0][0]
    dec = coords[1][0]
except IndexError:
    raise Warning('No candidates found.')
print(f'Hotspots:{coords}\n')
detection.closeXml()
tmodel = time.time() - t
print(f'Modelling: {tmodel} s\n')

analysis.usepnt = False
analysis.target = (ra, dec)
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
analysis.model = detname
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
tmodel += time.time() - t
print(f'Modelling2: {tmodel} s\n')

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
    raise Warning('No candidates found.')
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

logname = f'{rootpath}/DATA/outputs/crab/ctools3d_binned_blindfit.csv'
if first == 'True':
    hdr = 'texp sqrt_ts flux flux_err ra dec ttotal timport tsetup tobs tsky tblind tcube texpcube tpsfcube tbkgcube tmodel tfit tstat tflux\n'
    log = open(logname, 'w+')
    log.write(hdr)
    log.write(f'{texp} {np.sqrt(ts)} {phflux} {phflux_err} {ra} {dec} {ttotal} {timport} {tsetup} {tobs} {tsky} {tblind} {tcube} {texpcube} {tpsfcube} {tbkgcube} {tmodel} {tfit} {tstat} {tflux}\n')
    log.close()
else:
    log = open(logname, 'a')
    log.write(f'{texp} {np.sqrt(ts)} {phflux} {phflux_err} {ra} {dec} {ttotal} {timport} {tsetup} {tobs} {tsky} {tblind} {tcube} {texpcube} {tpsfcube} {tbkgcube} {tmodel} {tfit} {tstat} {tflux}\n')
    log.close()




