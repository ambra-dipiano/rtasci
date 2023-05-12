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
from rtasci.lib.RTACtoolsSimulation import make_obslist
from rtasci.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from rtasci.lib.RTAUtils import *
from rtasci.lib.RTAManageXml import ManageXml
from rtasci.lib.RTACtoolsBase import *
timport = time.time() - t
print(f'Imports : {timport} s\n')

t = time.time()
edisp = True
rootpath = "/data01/homes/cta/gammapy_integration"
obspath = f'{rootpath}/DATA/obs/crab/'
rtapath = f'{rootpath}/DATA/rta_products/crab/'
modelpath = f'{rootpath}/DATA/models/'
filename = f'{obspath}crab_offax_texp{texp}s_n01.fits'
obslist = filename.replace('.fits', '_obs.xml')
cube = filename.replace(obspath,rtapath).replace('.fits','_cube.xml')
expcube = cube.replace('.xml', '_exp.fits')
psfcube = cube.replace('.xml', '_psf.fits')
edispcube = cube.replace('.xml', '_edisp.fits')
bkgcube = cube.replace('.xml', '_bkg.fits')
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
analysis.caldb = 'prod3b-v2'
analysis.irf = 'South_z20_0.5h'
analysis.e = [0.05, 20]
analysis.t = [0, texp]
analysis.roi = 2.5
analysis.usepnt = True
analysis.input = obslist
analysis.output = cube
analysis.run_binning(prefix=cube.replace('.xml','_'), ebins=30, wbin=0.02)
tcube = time.time() - t
print(f'Binning: {tcube} s\n')

# exp cube
t = time.time()
analysis.input = obslist
analysis.output = expcube
if edisp:
    analysis.run_expcube(cube=cube.replace('.xml','_cta_01.fits'), ebins=30, wbin=0.02)
else:
    analysis.run_expcube(cube=cube.replace('.xml','_cta.fits'), ebins=30, wbin=0.02)
texpcube = time.time() - t
print(f'Exp. cube: {texpcube} s\n')

# psf cube
t = time.time()
analysis.input = obslist
analysis.output = psfcube
if edisp:
    analysis.run_psfcube(cube=cube.replace('.xml','_cta_01.fits'), ebins=30, wbin=0.5)
else:
    analysis.run_psfcube(cube=cube.replace('.xml','_cta.fits'), ebins=30, wbin=0.5)
tpsfcube = time.time() - t
print(f'Psf. cube: {tpsfcube} s\n')

if edisp:
    # edisp cube
    t = time.time()
    analysis.input = obslist
    analysis.output = edispcube
    analysis.run_edispcube(cube=cube.replace('.xml','_cta_01.fits'), ebins=30, wbin=0.5)
    tedispcube = time.time() - t
    print(f'Edisp. cube: {tedispcube} s\n')
else:
    tedispcube = np.nan
    edispcube=None

# bkg cube
t = time.time()
analysis.input = obslist
analysis.model = inmodel
analysis.output = bkgcube
if edisp:
    analysis.run_bkgcube(cube=cube.replace('.xml','_cta_01.fits'), model=outmodel)
else:
    analysis.run_bkgcube(cube=cube.replace('.xml','_cta.fits'), model=outmodel)
tbkgcube = time.time() - t
print(f'Bkg. cube: {tbkgcube} s\n')

# set model
t = time.time()
xml = ManageXml(outmodel)
xml.setTsTrue() 
xml.parametersFreeFixed(src_free=['Prefactor'], bkg_free=['Prefactor', 'Index'])
xml.closeXml()
tmodel = time.time() - t
print(f'Modelling: {tmodel} s\n')

# fitting
analysis.input = cube
analysis.model = outmodel
analysis.output = fitname
analysis.run_maxlikelihood(binned=True, edisp=edisp, exp=expcube, psf=psfcube, bkg=bkgcube, edispcube=edispcube)
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
print(f'pref={pref}')
print(analysis.e)
phflux = phflux_powerlaw(index, pref, pivot, analysis.e, unit='TeV')
phflux_err = phflux_powerlaw(index, err, pivot, analysis.e, unit='TeV')
print(f'PH-FLUX {phflux} +/- {phflux_err}\n')
tflux = time.time() - t
print(f'Flux points : {tflux} s\n')

ttotal = time.time() - clock0
print(f'Total time: {ttotal} s\n')
print('\n\n-----------------------------------------------------\n\n')

logname = f'{rootpath}/DATA/outputs/crab/ctools3d_binned_fit.csv'
row = f'{texp} {np.sqrt(ts)} {phflux} {phflux_err} {ttotal} {timport} {tsetup} {tobs} {tcube} {texpcube} {tpsfcube} {tedispcube} {tbkgcube} {tmodel} {tfit} {tstat} {tflux}\n'
if first == 'True':
    hdr = 'texp sqrt_ts flux flux_err ttotal timport tsetup tobs tcube texpcube tpsfcube tedispcube tbkgcube tmodel tfit tstat tflux\n'
    log = open(logname, 'w+')
    log.write(hdr)
    log.write(row)
    log.close()
else:
    log = open(logname, 'a')
    log.write(row)
    log.close()

print(row)

