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
from RTAscience.lib.RTAUtils import get_alert_pointing, get_mergermap
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
start_count = cfg.get('start_count')
trials = cfg.get('trials') + start_count
# sim parameters ---!
caldb = cfg.get('caldb')
irf = cfg.get('irf')
tobs = cfg.get('tobs')
onset = cfg.get('onset')
delay = cfg.get('delay')
tmax = tobs-onset+delay
emin = cfg.get('emin')
emax = cfg.get('emax')
roi = cfg.get('roi')
# conditions control ---!
set_ebl = cfg.get('set_ebl')

# paths ---!
datapath = cfg.get('data')

if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):  # obs parent folder
    os.mkdir(join(datapath, 'obs'))
# background model ---!
bkg_model = cfg.get('bkg')

# ---------------------------------------------------- TO-DO: filter population ---!!!

# ---------------------------------------------------- GRB population --- !!!
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
    pointing = get_alert_pointing(mergermap)
    for i in range(trials):
        count = start_count + i + 1
        name = f'ebl{count:06d}'
        # setup ---!
        sim = RTACtoolsSimulation()
        sim.configure(cfg)
        sim.seed = count
        sim.pointing = pointing
        sim.roi = roi
        #sim.e = [emin, emax]
        sim.tmax = tmax
        sim.caldb = caldb
        sim.irf = irf

        print(f'Simulate GRB + BKG with onset = {onset} s')
        sim.template = join(os.path.expandvars(cfg.get('catalog')).replace(cfg.get('data'), datapath), f'{runid}.fits')
        # load template ---!
        if not isfile(tcsv):
            raise FileExistsError(f'Table {csv} not found. Please set extrac_data=True')
        tbin_stop = sim.loadTemplate(source_name=runid, return_bin=True, data_path=join(datapath, f'extracted_data/{runid}'))

        event_bins = []
        # get time grid ---!
        sim.table = tcsv
        tgrid, tbin_start, tbin_stop = sim.getTimeSlices(GTI=(delay, tmax), return_bins=True) 
        # ----------------------------------------------- simulate ---!!!
        for i in range(tbin_stop-tbin_start-1):
            sim.t = [tgrid[i], tgrid[i + 1]]
            print(f'bin time interval = {sim.t}')
            sim.model = join(datapath, f'extracted_data/{runid}/{runid}_tbin{tbin_start+i:02d}.xml')
            event = join(grbpath, f'{name}_tbin{tbin_start+i:02d}.fits')
            event_bins.append(event)
            sim.output = event
            sim.run_simulation()

        # ---------------------------- merge in single photon list ---!!!
        phlist = join(grbpath, f'{name}.fits')
        sim.input = event_bins
        sim.output = phlist
        sim.appendEventsSinglePhList(GTI=[0, tobs])

        del sim
        print('remove bins')
        os.system('rm ' + join(grbpath, f'{name}*tbin*'))
