# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import argparse
from os.path import isdir, join, isfile
from rtasci.lib.RTAUtils import get_pointing
from rtasci.lib.RTAUtilsGW import get_offset
from rtasci.cfg.Config import Config
from rtasci.lib.RTAUtilsGW import get_alert_pointing_gw
from astropy import units as u
from astropy.coordinates import SkyCoord

parser = argparse.ArgumentParser(description='This script allows to find the offset between the target coordinate and the pointing coordinate.')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
args = parser.parse_args()
cfg = Config(args.cfgfile)

# GRB ---!
if cfg.get('runid') == 'all':
    runids = [f.replace('.fits', '') for f in os.listdir(cfg.get('catalog')) if isfile(join(cfg.get('catalog'), f))]
elif type(cfg.get('runid')) == str:
    runids = [cfg.get('runid')]
else:
    runids = cfg.get('runid')

# paths ---!
datapath = cfg.get('data')
if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')

# ------------------------------------------------------- loop runid --- !!!
for runid in runids:
    print(f'Processing runid: {runid}')
    # grb path ---!
    grbpath = join(datapath, 'obs', runid)  # folder that will host the phlist 
    if not isdir(grbpath):
        os.mkdir(grbpath)
    modelpath = join(datapath, f'extracted_data/{runid}')  # bin model folder
    if not isdir(modelpath):
        raise ValueError(f'Folder {runid} not found in {modelpath}')
    tcsv = join(datapath, f'extracted_data/{runid}/time_slices.csv')  # times table 
    if not isfile(tcsv):
        raise ValueError(f'Data from {runid} have not been correctly extracted.')
    mergerpath = os.path.expandvars(cfg.get('merger'))
    mergermap = get_mergermap(runid, mergerpath)
    template = f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits"

    # get alert pointing
    if type(cfg.get('offset')) == str and cfg.get('offset').lower() == 'gw':
        pointing = get_alert_pointing_gw(mergermap)
    else:
        pointing = list(get_pointing(template))
        if pointing[1] < 0:
            pointing[0] += 0.0
            pointing[1] += -cfg.get('offset')
        else:
            pointing[0] += 0.0
            pointing[1] += cfg.get('offset')

    print(f"Pointing: {pointing}")
    off = get_offset(template, mergermap)
    print(f"Offset: {off}")
    target = get_pointing(template)
    # find off-axis ---!
    offangle = SkyCoord(ra=target[0]*u.deg, dec=target[1]*u.deg, unit='deg', frame='icrs').separation(SkyCoord(ra=pointing[0]*u.deg, dec=pointing[1]*u.deg, unit='deg', frame='icrs')).deg
    print(f"Off-axis angle: {offangle}")