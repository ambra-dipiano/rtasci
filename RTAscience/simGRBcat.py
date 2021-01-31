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
from os.path import isdir, join, isfile
from lib.RTACtoolsSimulation import RTACtoolsSimulation
from lib.RTAManageXml import ManageXml
from lib.RTAUtils import get_alert_pointing

cfgfile = sys.argv[1]
pypath = str(os.path.dirname(os.path.abspath(__file__)))  
configuration = open(join(pypath, cfgfile) )
cfg = yaml.load(configuration, Loader=yaml.FullLoader)

# GRB ---!
if cfg['setup']['runid'] == 'all':
    runids = [f.replace('.fits', '') for f in os.listdir(cfg['path']['catalog']) if '_ebl' not in f and isfile(join(cfg['path']['catalog'], f))]
elif type(cfg['setup']['runid']) == str:
    runids = [cfg['setup']['runid']]
else:
    runids = cfg['setup']['runid']


# general ---!
trials = cfg['setup']['trials']
count = cfg['setup']['start_count']  
# sim parameters ---!
caldb = cfg['simulation']['caldb']
irf = cfg['simulation']['irf']
tobs = cfg['simulation']['tobs']
tonset = cfg['simulation']['tonset']
emin = cfg['simulation']['emin']
emax = cfg['simulation']['emax']
roi = cfg['simulation']['emax']
# conditions control ---!
set_ebl = cfg['options']['set_ebl']

# paths ---!
for path in cfg['path']:
    if '$' in path.values():
        datapath = os.path.expandvars(path.values())
datapath = cfg['path']['data']


breakpoint()
if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):  # obs parent folder
    os.mkdir(join(datapath, 'obs'))
# background model ---!
bkg_model = cfg['path']['bkg']

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
    merger = join(mergerpath.replace(cfg['path']['data'], datapath), 'f')

    # 
    pointing = get_alert_pointing()
    while count < trials:
        count += 1
        name = f'ebl{count:06d}'
        # setup ---!
        sim = RTACtoolsSimulation()
        sim.seed = count
        sim.nthreads = nthreads
        sim.pointing = pointing
        sim.roi = roi
        sim.e = [emin, emax]
        sim.tobs = tobs
        sim.caldb = caldb
        sim.irf = irf

        print(f'Simulate GRB + BKG with onset = {onset} s')
        sim.template = template
        sim.model = model_pl
        # add EBL to template ---!
        if set_ebl:
            print('Computing EBL absorption')
            sim.table = ebl_table  
            sim.zfetch = True
            sim.set_ebl = False
            if not isfile(template.replace('.fits', '_ebl.fits')):
                sim.addEBLtoFITS(template.replace('.fits', '_ebl.fits'), ext_name='EBL-ABS. SPECTRA')
            sim.set_ebl = set_ebl
            sim.template = template.replace('.fits', '_ebl.fits')
        # load template ---!
        if not isfile(tcsv):
            sim.extract_spectrum = True
        tbin_stop = sim.loadTemplate(source_name=runid, return_bin=True, data_path=join(datapath, f'extracted_data/{runid}'))

        event_bins = []
        # get time grid ---!
        sim.table = tcsv
        tgrid = sim.getTimeSlices(GTI=(0, tobs)) 
        # ----------------------------------------------- simulate ---!!!
        for i in range(tbin_stop):
            sim.t = [tgrid[i], tgrid[i + 1]]
            sim.model = join(datapath, f'extracted_data/{runid}/{runid}_tbin{i:02d}.xml')
            event = join(grbpath, f'{name}_tbin{i:02d}.fits')
            event_bins.append(event)
            sim.output = event
            sim.run_simulation()

        # ---------------------------- merge in single photon list ---!!!
        phlist = join(grbpath, f'{name}.fits')
        sim.input = event_bins
        sim.output = phlist
        sim.appendEventsSinglePhList(GTI=[0, tobs])

        del sim
        print('os.remove template bins')
        os.system('rm ' + join(grbpath, f'{name}*tbin*'))
