# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import sys
import astropy.units as u
import numpy as np
from os.path import join
from astropy.time import Time
from astropy.io import fits
from gammapy.data import DataStore

if len(sys.argv) < 2:
    raise ValueError('Missing nbins input from line command')
nbins = int(sys.argv[1])

path = os.path.join(os.path.expandvars('$DATA'), 'obs/crab/')
filename = 'crab_offax.fits'
nameroot = 'crab_offax'
fitsfile = join(path, filename)

data = DataStore.from_events_files([fitsfile])
observations = data.get_observations()
print(data.info())

# split 
start = observations[0].gti.time_start
duration = observations[0].gti.time_sum
times = start + np.linspace(0*u.s, duration, nbins)
time_intervals = [Time([tstart, tstop]) for tstart, tstop in zip(times[:-1], times[1:])]
bins = observations.select_time(time_intervals)
if not os.path.exists(path):
    os.mkdir(path)
for i, b in enumerate(bins):
    hdulist = fits.HDUList([fits.PrimaryHDU()])
    hdu = fits.BinTableHDU(b.events.table, name="EVENTS")
    hdulist.append(hdu)
    hdu = fits.BinTableHDU(b.gti.table, name="GTI")
    hdulist.append(hdu)
    hdulist.writeto(path / f"{nameroot}_texp{round(duration.value/nbins)}s_n{i+1:02d}.fits", overwrite=True)
    print(path / f"{nameroot}_texp{round(duration.value/nbins)}s_n{i+1:02d}.fits")