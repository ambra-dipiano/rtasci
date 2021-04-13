# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import sys
import argparse
import numpy as np
from time import time
from shutil import copy
from astropy.io import fits
from multiprocessing import Pool
from os.path import isdir, join, isfile
from RTAscience.cfg.Config import Config
from RTAscience.lib.RTAManageXml import ManageXml
from RTAscience.lib.RTACtoolsSimulation import RTACtoolsSimulation, make_obslist
from RTAscience.lib.RTAUtils import get_alert_pointing_gw, get_mergermap, get_pointing


def main(args):
    cfg = Config(args.cfgfile)
    # GRB ---!
    if cfg.get('runid') == 'all':
        runids = [f.replace('.fits', '') for f in os.listdir(cfg.get('catalog')) if isfile(join(cfg.get('catalog'), f))]
    elif type(cfg.get('runid')) == str:
        runids = [cfg.get('runid')]
    else:
        runids = cfg.get('runid')
    runids = sorted(runids)

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
        print(f"{'-'*50} #\nProcessing runid: {runid}")
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
        if mergermap == None:
            print(f'Skip runid {runid}. ')
            continue

        # get alert pointing
        if type(cfg.get('offset')) == str and cfg.get('offset').lower() == 'gw':
            pointing = get_alert_pointing_gw(mergermap)
        else:
            pointing = list(get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits"))
            if pointing[1] < 0:
                pointing[0] += 0.0
                pointing[1] += -cfg.get('offset')
            else:
                pointing[0] += 0.0
                pointing[1] += cfg.get('offset')

        # Dumping the Conf object to txt file
        dumpedConfig = os.path.join(grbpath, "config.yaml")
        if not os.path.isfile(dumpedConfig):
            copy(args.cfgfile, str(dumpedConfig))

        # ---------------------------------------------------- loop trials ---!!!
        if args.mp_enabled:                
            with Pool(args.mp_threads) as p:
                times = p.map(simulateTrial, [ (i, cfg, pointing, tmax, datapath, runid, tcsv, grbpath, bkg_model) for i in range(trials)])
        else:
            for i in range(trials):
                times = simulateTrial((i, cfg, pointing, tmax, datapath, runid, tcsv, grbpath, bkg_model))
        # time ---!
        if args.print.lower() == 'true':
            if len(times) > 1:
                print(f"Trial elapsed time (mean): {np.array(times).mean()}")
            else:
                print(f"Trial elapsed time: {times[0]}")    
        print('\n... done.\n')


def simulateTrial(trial_args):
    start_t = time()
    i=trial_args[0]
    cfg=trial_args[1]
    pointing=trial_args[2]
    tmax=trial_args[3]
    datapath=trial_args[4]
    runid=trial_args[5]
    tcsv=trial_args[6]
    grbpath=trial_args[7]
    bkg_model=trial_args[8]
    # initialise ---!
    count = cfg.get('start_count') + i + 1
    name = f'ebl{count:06d}'
    # setup ---!
    sim = RTACtoolsSimulation()
    sim.caldb = cfg.get('caldb')
    sim.irf = cfg.get('irf')
    sim.fov = cfg.get('roi')
    sim.e = [cfg.get('emin'), cfg.get('emax')]
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
    print(f'Simulate template seed={sim.seed}')
    for j in range(tbin_stop-tbin_start-1):
        sim.t = [tgrid[j]+cfg.get('onset'), tgrid[j + 1]+cfg.get('onset')]
        if args.print.lower() == 'true':
            print(f'GTI (bin) = {sim.t} s')
        sim.model = join(datapath, f'extracted_data/{runid}/{runid}_tbin{tbin_start+j:02d}.xml')
        event = join(grbpath, f'{name}_tbin{tbin_start+j:02d}.fits')
        event_bins.append(event)
        sim.output = event
        sim.run_simulation()
    # -------------------------------------------- shift time --- !!!
    if cfg.get('onset') != 0:
        if cfg.get('delay') != 0:
            raise ValueError('Bad configuration. Either "onset" or "delay" must be equal to 0.')
        # ------------------------------------ add background --- !!!
        print('Simulate bkg to append before the burst')
        bkg = os.path.join(grbpath, f'bkg{count:06d}.fits')
        event_bins.insert(0, bkg)
        sim.t = [0, cfg.get('onset')]
        if args.print.lower() == 'true':
            print(f"GTI (bkg) = {sim.t} s")
        sim.model = bkg_model
        sim.output = bkg
        sim.run_simulation()

    # ---------------------------------------- gather bins ---!!!
    if args.merge.lower() == 'true':
        print('Merge in photon-list')
        phlist = join(grbpath, f'{name}.fits')
        sim.input = event_bins
        sim.output = phlist
        sim.appendEventsSinglePhList(GTI=[cfg.get('delay'), cfg.get('delay')+cfg.get('tobs')])
        if args.print.lower() == 'true':
            h = fits.open(phlist)
            print('Check GTI and EVENTS time range:')
            print('************')
            print(h[2].data)
            print(h[1].data.field('TIME').min(), h[1].data.field('TIME').max())
            print('************')
            h.close()
    else:
        # observation list ---!
        obslist = join(grbpath, f'{name}.xml')
        if os.path.isfile(obslist):
            os.remove(obslist)
        make_obslist(obslist=obslist, items=event_bins, names=name)

    del sim
    if args.remove.lower() == 'true' and args.merge.lower() == 'true':
        # remove bins ---!
        os.system('rm ' + join(grbpath, f'{name}*tbin*'))
        if cfg.get('onset') != 0:
            # remove bkg ---!
            os.system('rm ' + join(grbpath, f'{name.replace("ebl", "bkg")}*'))
    # time ---!   
    elapsed_t = time()-start_t
    if args.print.lower() == 'true':
        print(f"Trial {count} took {elapsed_t} seconds")
    return (count, elapsed_t)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
    parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
    parser.add_argument('--merge', type=str, default='true', help='Merge in single phlist (true) or use observation library (false)')
    parser.add_argument('--remove', type=str, default='true', help='Keep only outputs')
    parser.add_argument('--print', type=str, default='false', help='Print out results')
    parser.add_argument('-mp', '--mp-enabled', type=str, default='false', help='To parallelize trials loop')
    parser.add_argument('-mpt', '--mp-threads', type=int, default=4, help='The size of the threads pool') 

    args = parser.parse_args()
    main(args)
