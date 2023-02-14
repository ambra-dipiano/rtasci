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
#from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import SkyCoord
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import EventList, GTI, Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel, PointSpatialModel
timport = time.time() - t
print(f'Imports : {timport} s\n')

t = time.time()
rootpath = str(os.path.dirname(os.path.abspath(__file__))).replace('cta-sag-sci/rtasci/timing', '')
caldb = f'{rootpath}/caldb/data/cta/prod3b-v2/bcf/South_z20_0.5h/irf_file.fits'
irfs = load_cta_irfs(caldb)
filename = f'{rootpath}/DATA/obs/crab/crab_offax_texp{texp}s_n01.fits'
obs_id = 1
print(f'Fits: {filename.replace(rootpath, "")}\n')
tsetup = time.time() - t
print(f'Setup : {tsetup} s\n')

# read phlist
t = time.time()
events = EventList.read(filename, hdu='EVENTS')
# get GTI
gti = GTI.read(filename, hdu='GTI')
# get pointing
pointing = events.pointing_radec
# create observation
observation = Observation.create(
    pointing=pointing, obs_id=f'{obs_id:02d}', tstart=gti.table['START'] * u.s, 
    tstop=gti.table['STOP'] * u.s, irfs=irfs, reference_time=gti.time_ref)
observation._events = events
#print(observation.gti)
observations = Observations() 
observations.append(observation)
# fix pointing info
observation.fixed_pointing_info
# target
#target = SkyCoord(pointing.ra, pointing.dec - 0.5 * u.deg, unit='deg', frame='icrs')
target = {'ra': 83.6331, 'dec': 22.0145}
tobs = time.time() - t
print(f'Create observation : {tobs} s\n')

# configure a 1d analysis
t = time.time()
config_1d = AnalysisConfig()
config_1d.general.log = {'level': 'warning'}
config_1d.datasets.type = "1d"
config_1d.datasets.stack = False
# define the ON region and make sure that PSF leakage is corrected
config_1d.datasets.on_region = dict(frame="icrs", lon='%s deg' %target['ra'], lat='%s deg' %target['dec'], radius='0.1 deg')
config_1d.datasets.containment_correction = True
# background
config_1d.datasets.background=dict(method="reflected", exclusion=None)
#config_1d.datasets.safe_mask.methods = ["edisp-bias"]
# define the energy binning for the spectra
config_1d.datasets.geom.axes.energy = dict(min='0.05 TeV', max='20 TeV', nbins=30)
config_1d.datasets.geom.axes.energy_true = dict(min='0.03 TeV', max='30 TeV', nbins=40)
config_1d.datasets.geom.selection.offset_max = '2.5 deg'
# fit
#config_1d.fit.fit_range = dict(min='0.03 TeV', max='30 TeV')
#config_1d.flux_points.energy = dict(min='0.03 TeV', max='30 TeV', nbins=3)
#config_1d.flux_points.source = 'Crab'
# write
#config_1d.write("config1d.yaml", overwrite=True)
#config_1d = AnalysisConfig.read("config1d.yaml")
tconf = time.time() - t
print(f'Configuration : {tconf} s\n')
print(target)
# instantiate data reduction passing directly the config object
t = time.time()
analysis_1d = Analysis(config_1d)
analysis_1d.observations = observations
analysis_1d.get_datasets()
tred = time.time() - t
print(f'Data Reduction : {tred} s\n')

# statistics
t = time.time()
stats = analysis_1d.datasets.info_table()
print(stats['sqrt_ts'], '\n')
tstat = time.time() - t
print(f'Statistics : {tstat} s\n')

# prepare models
t = time.time()
stacked_1d = analysis_1d.datasets.stack_reduce(name="stacked")
target = SkyCoord(target['ra'], target['dec'], unit='deg', frame='icrs')
spatial_model = PointSpatialModel(lon_0=target.ra, lat_0=target.dec, frame="icrs")
spectral_model = PowerLawSpectralModel(index=2.48, amplitude=2e-12 * u.Unit("1 / (cm2 s TeV)"))
spectral_model.parameters['index'].frozen = True
spatial_model.parameters['lon_0'].frozen = True
spatial_model.parameters['lat_0'].frozen = True
sky_model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model, name="Crab")
stacked_1d.models = sky_model
tmodel = time.time() - t
print(f'Modelling : {tmodel} s\n')

# fitting
t = time.time()
fit_1d = Fit([stacked_1d])
result_1d = fit_1d.run()
#print(result_1d.parameters.to_table(), '\n')
tfit = time.time() - t
print(f'Fitting : {tfit} s\n')

# flux
t = time.time()
phflux_err = spectral_model.integral_error(0.05 * u.TeV, 20 * u.TeV)
print(f'\nPH-FLUX {phflux_err.value[0]} +/- {phflux_err.value[1]}')
tflux = time.time() - t
print(f'\nFlux : {tflux} s\n')

ttotal = time.time() - clock0
print(f'Total time: {ttotal} s\n')
print('\n\n-----------------------------------------------------\n\n')

logname = f'{rootpath}/DATA/outputs/crab/gammapy1d_binned_fit.csv'
row = f'{texp} {stats["sqrt_ts"][0]} {phflux_err.value[0]} {phflux_err.value[1]} {ttotal} {timport} {tsetup} {tconf} {tred} {tstat} {tmodel} {tfit} {tflux}\n'
if first == 'True':
    hdr = 'texp sqrt_ts flux flux_err ttotal timport tsetup tobs tconf tred tstat tmodel tfit tflux\n'
    log = open(logname, 'w+')
    log.write(hdr)
    log.write(row)
    log.close()
else:
    log = open(logname, 'a')
    log.write(row)
    log.close()

print (row)