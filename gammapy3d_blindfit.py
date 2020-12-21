import time
t = time.time()
clock0 = time.time()
import astropy.units as u
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import EventList, GTI, Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapy.modeling import Fit
from gammapy.estimators import ExcessMapEstimator
from gammapy.estimators.utils import find_peaks
from gammapy.modeling.models import PointSpatialModel, PowerLawSpectralModel, SkyModel, FoVBackgroundModel
print(f'Imports : {time.time() - t} s\n')

t = time.time()
rootpath = '/home/ambra/Desktop/CTA/projects/'
caldb = f'{rootpath}/caldb/data/cta/prod3b-v2/bcf/South_z20_50h/irf_file.fits'
irfs = load_cta_irfs(caldb)
filename = f'{rootpath}/DATA/obs/crab/crab_onax.fits'
obs_id = 1
print(f'Setup : {time.time() - t} s\n')


# read phlist
t = time.time()
events = EventList.read(filename, hdu='EVENTS')
# get GTI
gti = GTI.read(filename, hdu='GTI')
# get pointing
pointing = events.pointing_radec
#print('Pointing :', pointing)
# create observation
observation = Observation.create(
    pointing=pointing, obs_id=f'{1:02d}', tstart=gti.table['START']*u.s, 
    tstop=gti.table['STOP']*u.s, irfs=irfs, reference_time=gti.time_ref)
observation._events = events
print(observation.gti)
observations = Observations()
observations.append(observation)
# fix pointing info
observation.fixed_pointing_info
print(f'Observation: {time.time() - t} s\n')

# configure a 3d analysis 
t = time.time()
config_3d = AnalysisConfig()
config_3d.general.log = {'level': 'warning'}
config_3d.observations.datastore = ''
config_3d.observations.obs_file = filename
# reduction type
config_3d.datasets.type = '3d'  # Analysis type is 3D
config_3d.datasets.stack = False  # We keep track of datasets in all bunches
# geometry of the map for 3d
config_3d.datasets.geom.wcs.skydir = {'lon': pointing.ra, 'lat': pointing.dec, 'frame': 'icrs'}  
config_3d.datasets.geom.wcs.fov = {'width': '6 deg', 'height': '6 deg'}
config_3d.datasets.geom.wcs.binsize = '0.02 deg'
# The FoV radius to use for cutouts
config_3d.datasets.geom.selection.offset_max = 2.5 * u.deg
# reconstructed energy axis for the counts map 
config_3d.datasets.geom.axes.energy = dict(min= "0.05 TeV", max="10 TeV", nbins=10)
# true energy axis for the IRF maps (should always be wider range and larger nbins)
config_3d.datasets.geom.axes.energy_true = dict(min= "0.03 TeV", max="30 TeV", nbins=30)
# backgroun
config_3d.datasets.background = {'method': 'fov_background', 'exclusion': None}
# safe mask from IRF and max offset
config_3d.datasets.safe_mask.methods = ['aeff-default', 'offset-max']
# what maps to compute
config_3d.datasets.map_selection = ['counts', 'exposure', 'background', 'psf', 'edisp']
# save the configuration for later and overwrite if already existing
#config_3d.write(filepath + 'tests/prototype3d.yaml', overwrite=True)
print(f'Configuration: {time.time() - t} s\n')

#print(config_3d)

# instantiate data reduction passing directly the config object
t = time.time()
analysis_3d = Analysis(config_3d)
# set observation (single - no list)
analysis_3d.observations = observations
# perform data reduction
analysis_3d.get_datasets()

print(analysis_3d.get_datasets())
print(f'Data Reduction {time.time() - t} s\n')

# stack and hotspot search
t = time.time()
stacked_3d = analysis_3d.datasets.stack_reduce(name="stacked_3d")
estimator = ExcessMapEstimator(correlation_radius='0.1 deg', selection_optional=[])
maps = estimator.run(stacked_3d)
hotspots_table = find_peaks(maps["sqrt_ts"].get_image_by_idx((0,)), threshold=5, min_distance='0.5 deg')
hotspots = SkyCoord(hotspots_table["ra"], hotspots_table["dec"])
print(f'Hotsposts: {hotspots}\n')
print(f'Blind search: {time.time() - t} s\n')

# target significance
t = time.time()
target = SkyCoord(hotspots.ra[0].deg, hotspots.dec[0].deg, unit='deg', frame='icrs')
target_region = CircleSkyRegion(target.icrs, 0.1 * u.deg)
stats = analysis_3d.datasets.info_table(cumulative=False)
print(stats['sqrt_ts'], '\n')
print(f'Statistics: {time.time() - t} s\n')

# modelling
t = time.time()
spatial_model = PointSpatialModel(lon_0=target.ra, lat_0=target.dec, frame="icrs")
spectral_model = PowerLawSpectralModel(index=2.3, amplitude=2e-12 * u.Unit("1 / (cm2 s TeV)"), reference=1 * u.TeV)
spectral_model.parameters['index'].frozen = True
spatial_model.parameters['lon_0'].frozen = True
spatial_model.parameters['lat_0'].frozen = True
sky_model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model, name="Crab")
bkg_model = FoVBackgroundModel(dataset_name="stacked_3d")
bkg_model.parameters['norm'].frozen = False
stacked_3d.models = [bkg_model, sky_model]
print(f'Modelling: {time.time() - t} s\n')

# fitting
t = time.time()
fit = Fit([stacked_3d])
result = fit.run()
print(result.parameters.to_table())
print(f'\nFitting : {time.time() - t} s\n')

# flux
t = time.time()
phflux = spectral_model.integral_error(0.05 * u.TeV, 20 * u.TeV)
enflux = spectral_model.energy_flux_error(0.05 * u.TeV, 20 * u.TeV)
print(f'\nPH-FLUX  {phflux[0]} +/- {phflux[1]}')
print(f'EN-FLUX {enflux[0]} +/- {phflux[1]}')

print(f'\nFlux points : {time.time() - t} s\n')
#print(result)

print(f'Total time: {time.time() - clock0} s\n')