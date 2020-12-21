import time
clock0 = time.time()
t = time.time()
#from astropy.coordinates import SkyCoord
from astropy import units as u
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import EventList, GTI, Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel
from gammapy.estimators import ExcessMapEstimator
from gammapy.estimators.utils import find_peaks
print(f'Imports : {time.time() - t} s\n')

t = time.time()
rootpath = '/home/ambra/Desktop/CTA/projects/REMOTE/'
caldb = f'{rootpath}/caldb/data/cta/prod3b-v2/bcf/South_z20_50h/irf_file.fits'
irfs = load_cta_irfs(caldb)
filename = f'{rootpath}/gammapy_integration/DATA/crab/crab_texp10s_n00.fits'
obs_id = 1
print(f'Setup : {time.time() - t} s\n')

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
print(observation.gti)
observations = Observations() 
observations.append(observation)

# fix pointing info
observation.fixed_pointing_info
# target
#target = SkyCoord(pointing.ra, pointing.dec - 0.5 * u.deg, unit='deg', frame='icrs')
target = {'ra': pointing.ra.value, 'dec': pointing.dec.value - 0.5}
print(f'Create observation : {time.time() - t} s\n')

# configure a 1d analysis ---- REQUIRES 3D + 1D schould try mixing ctools and gammapy 
t = time.time()
config_1d = AnalysisConfig.read("config1d.yaml")
print(f'Configuration : {time.time() - t} s\n')

#print(config_1d)

# instantiate data reduction passing directly the config object
t = time.time()
analysis_1d = Analysis(config_1d)
analysis_1d.observations = observations
analysis_1d.get_datasets()
print(f'Data Reduction : {time.time() - t} s\n')

# statistics
t = time.time()
stats = analysis_1d.datasets.info_table()
print(stats['sqrt_ts'])
print(f'Statistics : {time.time() - t} s\n')

# stack and hotspot search
t = time.time()
stacked_1d = analysis_1d.datasets.stack_reduce(name="stacked")
estimator = ExcessMapEstimator(correlation_radius='0.1 deg', selection_optional=[])
maps = estimator.run(stacked_1d)
hotspots_table = find_peaks(maps["sqrt_ts"].get_image_by_idx((0,)), threshold=5, min_distance='0.5 deg')
hotspots = SkyCoord(hotspots_table["ra"], hotspots_table["dec"])
print('hotspot', hotspots, '\npointing', pointing, '\ntarget', target, '\ntarget region', target_region)
print(f'Blind search: {time.time() - t} s\n')

# modelling
t = time.time()
spectral_model = PowerLawSpectralModel(index=2.3, amplitude="2e-12 TeV-1 cm-2 s-1")
spectral_model.parameters['index'].frozen = True
sky_model=SkyModel(spectral_model=spectral_model, name="Crab")
stacked_1d.models = sky_model
print(f'Modelling : {time.time() - t} s\n')

# hot-spots
estimator = ExcessMapEstimator(correlation_radius='0.1 deg', selection_optional=[])
maps = estimator.run(stacked_1d)
hotspots_table = find_peaks(maps["sqrt_ts"].get_image_by_idx((0,)), threshold=5, min_distance='0.1 deg')
hotspots = SkyCoord(hotspots_table["ra"], hotspots_table["dec"])
print('TS MAX', maps["sqrt_ts"].max())
print('COORDS', hotspots.ra.deg, hotspots.dec.deg)






# fitting
t = time.time()
fit_1d = Fit([stacked_1d])
result_1d = fit_1d.run()
#print(result_1d)
print(result_1d.parameters.to_table())
print(f'Fitting : {time.time() - t} s\n')

# flux
t = time.time()
print(f'PH-FLUX {spectral_model.integral(0.05 * u.TeV, 20 * u.TeV)} +/- {spectral_model.integral_error(0.05 * u.TeV, 20 * u.TeV)}')
print(f'EN-FLUX {spectral_model.energy_flux(0.05 * u.TeV, 20 * u.TeV)} +/- {spectral_model.energy_flux_error(0.05 * u.TeV, 20 * u.TeV)}')

print(f'Flux points : {time.time() - t} s\n')

print(f'Total : {time.time() - clock0} s\n')