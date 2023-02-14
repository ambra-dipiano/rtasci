# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import numpy as np
import healpy as hp
from rtasci.lib.RTAUtils import get_pointing

# retrieve telescope pointing coordinates from alert probability map ---!
def get_alert_pointing_compressed(merger_map):
    '''Given a compressed merger map, it returns the maximum sky localisation probability coordinates.'''
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


# retrieve telescope pointing coordinates from alert probability map ---!
def get_alert_pointing_gw(merger_map):
    '''Given a merger map, it returns the maximum sky localisation probability coordinates.'''
    if not os.path.isfile(merger_map):
        raise ValueError(f'Merger map {merger_map} not found.')
    # load map ---!
    map = hp.read_map(merger_map, dtype=None)
    pixels = len(map)
    axis = hp.npix2nside(pixels)
    # search max prob coords ---!
    pmax = np.argmax(map)
    theta, phi = hp.pix2ang(axis, pmax)
    pointing = (np.rad2deg(phi), np.rad2deg(0.5 * np.pi - theta))
    return pointing

# center of fov from FITS ---!
def get_offset(fits_file, merger_map):
    '''Given a template and a merger map, returns the offset (in RA and DEC) between the maximum sky localisation probability and the target coordinates.'''
    true = get_pointing(fits_file)
    alert = get_alert_pointing_gw(merger_map)
    return (true[0]-alert[0], true[1]-alert[1])