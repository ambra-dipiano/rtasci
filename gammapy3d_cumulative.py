import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from pathlib import Path
from astropy.time import Time
from regions import CircleSkyRegion
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import EventList, GTI, Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapy.modeling import Fit
from gammapy.estimators import ExcessMapEstimator
from gammapy.estimators.utils import find_peaks
from gammapy.modeling.models import (
    PointSpatialModel,
    PowerLawSpectralModel,
    SkyModel,
    FoVBackgroundModel,
    Models,
)

caldb = "/data01/homes/cta/caldb/data/cta/prod3b-v2/bcf/South_z20_50h/irf_file.fits"
irfs = load_cta_irfs(caldb)

path = Path("/data01/homes/cta/gammapy_integration/tests/TEST_DATA/bins")
obs_id = 1

# collect bins in observation
observations = Observations() 
for i in range(10):
    events = EventList.read(path / f"data_{obs_id}_bunch_{i:02d}.fits", hdu="EVENTS")
    gti = GTI.read(path / f"data_{obs_id}_bunch_{i:02d}.fits", hdu="GTI")
    pointing = events.pointing_radec
    observation = Observation.create(
        pointing=pointing,
        obs_id=f'{i+1:02d}',
        tstart=gti.table["START"]*u.s,
        tstop=gti.table["STOP"]*u.s,
        irfs=irfs,
        reference_time=gti.time_ref
    )
    observation._events = events
    observations.append(observation)

print('Observations IDs :', observations.ids)

# get pointing (same for all bins)
pointing = observation.pointing_radec

print('Pointing :', pointing)

# configure a 3d analysis