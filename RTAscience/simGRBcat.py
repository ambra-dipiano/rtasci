# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import numpy as np
import os
import sys
import argparse
from os.path import isdir, join, isfile
from RTAscience.lib.RTACtoolsSimulation import RTACtoolsSimulation
from RTAscience.lib.RTAManageXml import ManageXml
from RTAscience.lib.RTAUtils import get_alert_pointing, get_mergermap, get_pointing
from RTAscience.cfg.Config import Config

parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
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

# general ---!
trials = cfg.get('trials')
tmax = cfg.get('tobs')-cfg.get('onset')+cfg.get('delay')

# paths ---!
datapath = cfg.get('data')
if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):  # obs parent folder
    os.mkdir(join(datapath, 'obs'))
# background model ---!
bkg_model = cfg.get('bkg')

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

    # get alert pointing
    if type(cfg.get('offset')) == str and cfg.get('offset').lower() == 'gw':
        pointing = get_alert_pointing(mergermap)
    else:
        pointing = list(get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits"))
        if pointing[1] < 0:
            pointing[0] += 0.0
            pointing[1] += -cfg.get('offset')
        else:
            pointing[0] += 0.0
            pointing[1] += cfg.get('offset')

    # ---------------------------------------------------- loop trials ---!!!
    for i in range(trials):
        count = cfg.get('start_count') + i + 1
        name = f'ebl{count:06d}'
        # setup ---!
        sim = RTACtoolsSimulation()
        sim.configure(cfg)
        sim.seed = count
        sim.set_ebl = cfg.get('set_ebl')
        sim.pointing = pointing
        sim.tmax = tmax

        # get time grid ---!
        sim.template = join(os.path.expandvars(cfg.get('catalog')).replace(cfg.get('data'), datapath), f'{runid}.fits')
        event_bins = []
        sim.table = tcsv
        tgrid, tbin_start, tbin_stop = sim.getTimeSlices(GTI=(cfg.get('delay'), tmax), return_bins=True) 

        # -------------------------------------------------------- simulate ---!!!
        for j in range(tbin_stop-tbin_start-1):
            sim.t = [tgrid[j], tgrid[j + 1]]
            print(f'GTI (bin) = {sim.t}')
            sim.model = join(datapath, f'extracted_data/{runid}/{runid}_tbin{tbin_start+j:02d}.xml')
            event = join(grbpath, f'{name}_tbin{tbin_start+j:02d}.fits')
            event_bins.append(event)
            sim.output = event
            sim.run_simulation()

        # ---------------------------------------- merge in single photon list ---!!!
        print('Merge bins in photon-list')
        phlist = join(grbpath, f'{name}.fits')
        sim.input = event_bins
        sim.output = phlist
        sim.appendEventsSinglePhList(GTI=[0, cfg.get('tobs')])

        del sim
        print('Remove bins')
        os.system('rm ' + join(grbpath, f'{name}*tbin*'))
print('... done.\n')

