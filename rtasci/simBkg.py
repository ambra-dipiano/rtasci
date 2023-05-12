# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import numpy as np
import os, argparse
import time
from time import time
from shutil import copy
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool
from os.path import isdir, expandvars
from rtasci.cfg.Config import Config
from rtasci.lib.RTACtoolsSimulation import RTACtoolsSimulation
from rtasci.lib.RTAUtils import get_pointing, get_mergermap, str2bool
from rtasci.lib.RTAUtilsGW import get_alert_pointing_gw

def main(args):
    cfg = Config(args.cfgfile)
    # general ---!
    if cfg.get('simtype') != 'bkg':
        raise ValueError('This script only allows bakground simulations')
    trials = cfg.get('trials')  # trials
    tobs = cfg.get('tobs')  # total obs time (s)
    runid = cfg.get('runid')
    if type(runid) != str or runid == 'all':
        raise ValueError('Currently only 1 runid is allowed.')

    # paths ---!
    datapath = cfg.get('data')
    if not isdir(datapath):  # main data folder
        raise ValueError('Please specify a valid path')
    if args.output_dir:
        bkgpath = Path(args.output_dir)
    else:
        bkgpath = Path(datapath).joinpath('obs', 'backgrounds')
    bkgpath.mkdir(parents=True, exist_ok=True)

    # background model ---!
    bkg_model = expandvars(cfg.get('bkg'))  # XML background model

    # ------------------------------------------------------- loop runid --- !!!
    # get alert pointing
    if type(cfg.get('offset')) == str and cfg.get('offset') == 'gw':
        mergerpath = os.path.expandvars(cfg.get('merger'))
        mergermap = get_mergermap(runid, mergerpath)
        if mergermap == None:
            raise ValueError(f'Merger map of runid {runid} not found. ')
        pointing = get_alert_pointing_gw(mergermap)
    else:
        if runid == 'crab':
            pointing = [83.6331, 22.0145]
        else:
            pointing = list(get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits"))
        if pointing[1] < 0:
            pointing[0] += 0.0
            pointing[1] += -cfg.get('offset')
        else:
            pointing[0] += 0.0
            pointing[1] += cfg.get('offset')

    # Dumping the Conf object to txt file
    dumpedConfig = os.path.join(bkgpath, "config.yaml")
    if not os.path.isfile(dumpedConfig):
        copy(args.cfgfile, str(dumpedConfig))

    # ---------------------------------------------------- loop trials ---!!!
    if args.mp_enabled:
        with Pool(args.mp_threads) as p:
            times = p.map(simulateTrial, [ (i, cfg, pointing, bkg_model, bkgpath, tobs, args.remove) for i in range(trials)])
    else:
        for i in range(trials):
            times = simulateTrial((i, cfg, pointing, bkg_model, bkgpath, tobs, args.remove))
    # time ---!
    if args.print:
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
    bkg_model=trial_args[3]
    bkgpath=trial_args[4]
    tobs=trial_args[5]
    remove_logs=trial_args[6]
    # initialise ---!
    count = cfg.get('start_count') + i + 1
    name = f'bkg{count:08d}'
    # setup ---!
    sim = RTACtoolsSimulation()
    sim.seed = count
    sim.pointing = pointing
    sim.caldb = cfg.get('caldb')
    sim.irf = cfg.get('irf')
    sim.fov = cfg.get('roi')
    sim.e = [cfg.get('emin'), cfg.get('emax')]


    print(f"[{datetime.now().strftime('%d/%m/%Y_%H:%M:%S')}] Simulate empty fields for runid = {cfg.get('runid')} with seed = {count}", flush=True)
    sim.seed = count
    sim.t = [0, tobs]
    bkg = os.path.join(bkgpath, f'{name}.fits')
    sim.model = bkg_model
    sim.output = bkg
    if args.print:
        print(f"Simulation {bkg}")    
    sim.run_simulation()
    if remove_logs:
        Path(sim.output).with_suffix('.log').unlink()
    sim.input = bkg
    sim.sortObsEvents()
    del sim
    # time ---!
    elapsed_t = time()-start_t
    if args.print:
        print(f"Trial {count} took {elapsed_t} seconds.")
    print(f".. done [{datetime.now().strftime('%d/%m/%Y_%H:%M:%S')}]", flush=True)
    return (count, elapsed_t)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Simulate empty fields.')
    parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
    parser.add_argument('--print', type=str2bool, default='false', help='Print out results')
    parser.add_argument('--remove', type=str2bool, default='true', help='Keep only .fits files and not .log')
    parser.add_argument('-mp', '--mp-enabled', type=str2bool, default='false', help='To parallelize trials loop')
    parser.add_argument('-mpt', '--mp-threads', type=int, default=4, help='The size of the threads pool')
    parser.add_argument('-out', '--output-dir', type=str, required=False, default="", help='The path to the output directory')
    args = parser.parse_args()

    main(args)
