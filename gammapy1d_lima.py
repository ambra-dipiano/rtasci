import time
clock0 = time.time()
t = time.time()
#from astropy.coordinates import SkyCoord
from astropy import units as u
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import EventList, GTI, Observation, Observations
from gammapy.irf import load_cta_irfs
print(f'Imports : {time.time() - t} s')

t = time.time()
rootpath = '/home/ambra/Desktop/CTA/projects/REMOTE/'
caldb = f'{rootpath}/caldb/data/cta/prod3b-v2/bcf/South_z20_50h/irf_file.fits'
irfs = load_cta_irfs(caldb)
filename = f'{rootpath}/gammapy_integration/DATA/crab/crab_texp2s_n00.fits'
obs_id = 1
print(f'Setup : {time.time() - t} s')

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
print('\n', observation.gti)
observations = Observations() 
observations.append(observation)

# fix pointing info
observation.fixed_pointing_info
# target
#target = SkyCoord(pointing.ra, pointing.dec - 0.5 * u.deg, unit='deg', frame='icrs')
target = {'ra': pointing.ra.value, 'dec': pointing.dec.value - 0.5}
print(f'Create observation : {time.time() - t} s')

# configure a 1d analysis
t = time.time()
config_1d = AnalysisConfig()
config_1d.general.log = {'level': 'warning'}
config_1d.datasets.type = "1d"
config_1d.datasets.stack = False
# define the ON region and make sure that PSF leakage is corrected
config_1d.datasets.on_region = dict(frame="icrs", lon='%s deg' %target['ra'], lat='%s deg' %target['dec'], radius='0.2 deg')
config_1d.datasets.containment_correction = True
# background
config_1d.datasets.background=dict(method="reflected", exclusion=None)
config_1d.datasets.safe_mask.methods = ["edisp-bias"]
#config_1d.datasets.safe_mask.parameters = dict(edisp_bias=10)
# define the energy binning for the spectra
config_1d.datasets.geom.axes.energy = dict(min='0.05 TeV', max='20 TeV', nbins=20)
config_1d.datasets.geom.axes.energy_true = dict(min='0.03 TeV', max='30 TeV', nbins=30)
config_1d.datasets.geom.selection.offset_max = '3 deg'
# write
#config_1d.write(f"{rootpath}/cfg/config1d.yaml", overwrite=True)
#config_1d = AnalysisConfig.read(f"{rootpath}/cfg/config1d.yaml")
print(f'Configuration : {time.time() - t} s')

#print(config_1d)

# instantiate data reduction passing directly the config object
t = time.time()
analysis_1d = Analysis(config_1d)
analysis_1d.observations = observations
analysis_1d.get_datasets()
print(f'Data Reduction : {time.time() - t} s')

# statistics
t = time.time()
stats = analysis_1d.datasets.info_table()
#print(stats.columns)
print('\n', stats['sqrt_ts'])
print('\n', stats['excess'])
print(f'Statistics : {time.time() - t} s')




print(f'Total : {time.time() - clock0} s')