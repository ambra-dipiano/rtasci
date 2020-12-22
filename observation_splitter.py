from pathlib import Path
import astropy.units as u
import numpy as np
from astropy.time import Time
import astropy.io.fits as fits
from gammapy.data import DataStore
import os

pathin = Path('/home/ambra/Desktop/CTA/projects/DATA/selections/crab/')
pathout = Path('/home/ambra/Desktop/CTA/projects/DATA/selections/crab')
filename = 'crab_offax_texp10s_n01.fits'
nameroot = 'crab_offax'
fitsfile = pathin / filename
nbins = 5

data_store = DataStore.from_events_files([fitsfile])
observations = data_store.get_observations()
print(data_store.info())

# split 
start = observations[0].gti.time_start
duration = observations[0].gti.time_sum
times = start + np.linspace(0*u.s, duration, nbins)
time_intervals = [Time([tstart, tstop]) for tstart, tstop in zip(times[:-1], times[1:])]
bunches = observations.select_time(time_intervals)
if not os.path.exists(pathout):
    os.mkdir(pathout)
pathout.mkdir(exist_ok=True)
for i, bunch in enumerate(bunches):
    hdulist = fits.HDUList([fits.PrimaryHDU()])
    hdu = fits.BinTableHDU(bunch.events.table, name="EVENTS")
    hdulist.append(hdu)
    hdu = fits.BinTableHDU(bunch.gti.table, name="GTI")
    hdulist.append(hdu)
    hdulist.writeto(pathout / f"{nameroot}_texp{round(duration.value/nbins)}s_n{i+1:02d}.fits", overwrite=True)
    print(pathout / f"{nameroot}_texp{round(duration.value/nbins)}s_n{i+1:02d}.fits")