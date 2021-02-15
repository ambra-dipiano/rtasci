# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import astropy.units as u
import healpy as hp
import numpy as np
from astropy.io import fits

# center of fov from FITS ---!
def get_pointing(fits_file):
    with fits.open(fits_file) as hdul:
        ra = hdul[0].header['RA']
        dec = hdul[0].header['DEC']
    return (ra, dec)

# retrieve telescope pointing coordinates from alert probability map ---!
def get_alert_pointing(merger_map):
    # load map ---!
    merger_map = merger_map.replace('.gz','')
    os.system(f'gunzip {merger_map}.gz')
    map = hp.read_map(merger_map, dtype=None)
    pixels = len(map)
    axis = hp.npix2nside(pixels)
    # search max prob coords ---!
    pmax = np.argmax(map)
    theta, phi = hp.pix2ang(axis, pmax)
    pointing = (np.rad2deg(phi), np.rad2deg(0.5 * np.pi - theta))
    os.system(f'gzip {merger_map}')
    return pointing

def increase_exposure(x, function='double'):
    y = None
    if function.lower() == 'double':
        y = x*2
    elif function.lower() == 'power2':
        y = x**2
    elif function.lower() == 'times10':
        y = x*10
    return y

# compute integral photon flux for PL model ---!
def phflux_powerlaw(gamma, k0, e0=1, erange=(0.03, 150.0), unit='TeV'):
    if unit == 'eV':
        conv = 1e-6
    elif unit == 'keV':
        conv = 1e-3
    elif unit == 'MeV':
        conv = 1
    elif unit == 'GeV':
        conv = 1e3
    else:
        conv = 1e6
    e1 = erange[0] * conv
    e2 = erange[1] * conv
    delta = gamma + 1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

# compute integral energy flux for PL model ---!
def enflux_powerlaw(gamma, k0, e0=1, erange=(0.03, 150.0), unit='TeV'):
    if unit == 'eV':
        conv = 1.60218e-12
    elif unit == 'keV':
        conv = 1.60218e-9
    elif unit == 'MeV':
        conv = 1.60218e-6
    elif unit == 'GeV':
        conv = 1.60218e-3
    else:
        conv = 1.60218
    e1 = erange[0] * conv
    e2 = erange[1] * conv
    k0 *= conv
    e0 *= conv
    delta = gamma+1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

# returns a random total delay time (slew time + gw latency) within given ranges ---!
def totalDelay(slew=(0,50), gw_latency=(0,36000)):
    tslew = np.random.uniform(slew[0], slew[1], 1)
    tgw = np.random.uniform(gw_latency[0], gw_latency[1])
    delay = tslew + tgw
    return delay

# get wobble pointing for given target
def wobble_pointing(target, nrun, clockwise=True, offset=0.5):
    if clockwise:
        wobble = [(0., offset), (-offset, 0.), (0., -offset), (offset, 0.)]
    else:
        wobble = [(0., offset), (offset, 0.), (0., -offset), (-offset, 0.)] 
    wobble_index = nrun % 4
    print(wobble[wobble_index])
    pointing = (target[0] + wobble[wobble_index][0], target[1] + wobble[wobble_index][1])
    return pointing

# get mergermap file
def get_mergermap(run, path):
    runid = run.split('_')
    merger = f'{runid[0]}_Merger{runid[1]}_skymap.fits.gz'
    if merger in os.listdir(path):
        return os.path.join(path, merger)
    else:
        raise FileExistsError(f'Merger map {merger} for not found in path: {path}')
