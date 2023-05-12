# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <.dipiano@inaf.it>
# *******************************************************************************

from sagsci.tools.utils import *
from sagsci.wrappers.rtaph.photometry import Photometrics, aeff_eval

# some parameters
pointing = {'ra': 57, 'dec': 33}
region_radius = 0.2
offset = 0.5
max_offset = 2
livetime = 10
args = 'the other required stuff'

# example of data collection
data = {'region': [], 'counts':[], 'flux':[]}

# define the rings in FOV
while offset <= max_offset:
    # define the starting (aka target not necessary source) region
    target_region = {'ra': pointing[0]+offset, 'dec': pointing[1], 'rad': region_radius}

    # init photometry
    phm = Photometrics(*args)
    # get response for current offset
    region_eff_resp = aeff_eval(target_region, pointing)
    # get ring of regions from starting position
    ring_regions = phm.find_off_regions(target_region, pointing, method='reflection', *args)
    # add starting position region to dictionary 
    ring_regions += target_region

    # for each region in ring count events and get flux
    for region in ring_regions:
        # counts in region within time window
        counts = phm.region_counter(region, *args)
        # normalised value (flux) in region within time window
        flux = counts / livetime / region_eff_resp 
        # save data somewhere
        data['region'].append(region)
        data['counts'].append(counts)
        data['flux'].append(flux)

    # increment offset to get next offset ring
    offset += region_radius*2
