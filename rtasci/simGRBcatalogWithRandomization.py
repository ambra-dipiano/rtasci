# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
# *******************************************************************************

import os
import random
import argparse
import numpy as np
from time import time
from pathlib import Path
from astropy.io import fits
from datetime import datetime
from multiprocessing import Pool
from os.path import join, isfile
from shutil import move, rmtree
from rtasci.cfg.Config import Config
from rtasci.lib.RTACtoolsSimulation import RTACtoolsSimulation, make_obslist
from rtasci.lib.RTAUtils import get_alert_pointing_gw, get_mergermap, get_pointing, str2bool

class TrialOutput:
    def __init__(self, trial_id, run_id, elapsed_time, errors=False, error_msg=''):
        self.trial_id = trial_id
        self.run_id = run_id
        self.elapsed_time = elapsed_time
        self.errors = errors
        self.error_msg = error_msg
    def __str__(self):
        return f"Trial {self.trial_id} (run {self.run_id}) took {self.elapsed_time:.2f} seconds. Error message: {self.error_msg}"

def is_randomizable(val):
    if isinstance(val, str) and "random" in val:
        return True
    return False

def randomize(val):
    _min, _max = val.replace("random$","").split("-")
    return round(random.uniform(float(_min), float(_max)), 2)

def config_runid(cfg):
    # GRB ---!
    if cfg.get('runid') == 'all':
        runids = [f.replace('.fits', '') for f in os.listdir(cfg.get('catalog')) if isfile(join(cfg.get('catalog'), f))]
    elif type(cfg.get('runid')) == str:
        runids = [cfg.get('runid')]
    else:
        runids = cfg.get('runid')
    return sorted(runids)

def create_output_dirs(output_dir, cfg, runids):
    start = time()
    # Create output dirs and dump the configuration file inside
    for runid in runids:
        output_path = Path(output_dir).joinpath(runid)
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        cfg.set("runid", runid) # it could be a list or "all" but I want a single runid
        cfg.dump(output_path.joinpath("config.yaml")) 
    print(f"Output dirs created in {round(time()-start, 3)} seconds")

def now():
    return datetime.now().strftime('%d/%m/%Y_%H:%M:%S')

def get_pointing_utility(runid, cfg):

    # get alert pointing
    pointing = None
    offset = cfg.get('offset')

    if isinstance(offset, str) and offset.lower() == 'gw':
        mergerpath = os.path.expandvars(cfg.get('merger'))
        mergermap = get_mergermap(runid, mergerpath)
        pointing = get_alert_pointing_gw(mergermap)
    
    else:
        pointing = list(get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits"))
        if pointing[1] < 0:
            pointing[0] += 0.0
            pointing[1] += -offset
        else:
            pointing[0] += 0.0
            pointing[1] += offset

    return pointing

def simulate_trial_bkg(input_args):
    
    runid, trial_id, cfg, args = input_args

    start_t = time()

    sim_output_path = args.output_dir.joinpath("backgrounds") 
    sim_output_path.mkdir(parents=True, exist_ok=True)

    offset = cfg.get('offset')
    tobs = cfg.get('tobs')
    bkg_model = cfg.get('bkg')

    if args.print:
        print(f"start_count: {cfg.get('start_count')} trial id: {trial_id}")
    
    count = cfg.get('start_count') + trial_id + 1
    name = f'runid_notemplate_trial_{count:010d}_simtype_{cfg.get("simtype")}_onset_0_delay_0_offset_{offset}_tobs_{tobs}.fits'

    try:
        pointing = get_pointing_utility(runid, cfg)
    except Exception as e:
        return TrialOutput(trial_id, runid, time()-start_t, True, str(e))

    # setup ---!
    sim = RTACtoolsSimulation()
    sim.seed = count
    sim.pointing = pointing
    sim.caldb = cfg.get('caldb')
    sim.irf = cfg.get('irf')
    sim.fov = cfg.get('roi')
    sim.e = [cfg.get('emin'), cfg.get('emax')]

    if args.print:
        print(f"[{now()}] Simulate empty fields for runid = {cfg.get('runid')} with seed = {count}")
    
    sim.seed = count
    sim.t = [0, tobs]
    bkg = str(sim_output_path.joinpath(name))
    sim.model = bkg_model
    sim.output = bkg
    sim.run_simulation()
    #if remove_logs:
    #    Path(sim.output).with_suffix('.log').unlink()
    sim.input = bkg
    sim.sortObsEvents()
    del sim
    # time ---!
    elapsed_t = time()-start_t
    if args.print:
        print(f"Trial {count} took {elapsed_t} seconds.")
        print(f".. done [{now()}]", flush=True)

    # time ---!   
    elapsed_t = time()-start_t

    return TrialOutput(trial_id, runid, elapsed_t)

def simulate_trial_grb(input_args):

    runid, trial_id, cfg, args = input_args

    start_t = time()

    sim_output_path_base = args.output_dir.joinpath(runid) 
    sim_output_path = sim_output_path_base.joinpath(f"trial_{trial_id}")  # folder that will host the phlist 
    
    if not sim_output_path.exists():
        sim_output_path.mkdir(parents=True, exist_ok=True)

    datapath = Path(cfg.get('data'))

    modelpath = datapath.joinpath("extracted_data", runid) # bin model folder
    if not modelpath.exists():
        return TrialOutput(trial_id, runid, 0, True, f'Folder {runid} not found in {modelpath}')

    tcsv = datapath.joinpath("extracted_data", runid, "time_slices.csv") # times table 
    if not tcsv.exists():
        return TrialOutput(trial_id, runid, time()-start_t, True, f"File {tcsv} not found")

    onset  = cfg.get('onset')
    delay  = cfg.get('delay')
    offset = cfg.get('offset')

    # add more if needed
    if(is_randomizable(cfg.get('onset'))):
        onset = randomize(cfg.get('onset'))
        
    if(is_randomizable(cfg.get('delay'))):
        delay = randomize(cfg.get('delay'))

    if(is_randomizable(cfg.get('offset'))):
        offset = randomize(cfg.get('offset'))

    try:
        pointing = get_pointing_utility(runid, cfg)
    except Exception as e:
        return TrialOutput(trial_id, runid, time()-start_t, True, str(e))

    tmax = cfg.get('tobs')- onset + delay
    bkg_model = cfg.get('bkg')

    # initialise ---!
    print(f"start_count: {cfg.get('start_count')} trial id: {trial_id}")
    count = cfg.get('start_count') + trial_id + 1
    name = f'runid_{runid}_trial_{count:010d}_simtype_{cfg.get("simtype")}_onset_{onset}_delay_{delay}_offset_{offset}'
    
    # setup ---!
    sim = RTACtoolsSimulation()
    if type(cfg.get('caldb')) == list:
        sim.caldb = cfg.get('caldb')[0]
    else:
        sim.caldb = cfg.get('caldb')
    if type(cfg.get('irf')) == list:
        sim.irf = cfg.get('irf')[0]
    else:
        sim.irf = cfg.get('irf')
    sim.fov = cfg.get('roi')
    sim.e = [cfg.get('emin'), cfg.get('emax')]
    sim.seed = count
    sim.set_ebl = cfg.get('set_ebl')
    sim.pointing = pointing
    if args.print:
        print(f'Pointing = {sim.pointing} s')
    sim.tmax = tmax

    # get time grid ---!
    sim.template = join(os.path.expandvars(cfg.get('catalog')).replace(cfg.get('data'), str(datapath)), f'{runid}.fits')
    event_bins = []
    sim.table = tcsv
    tgrid, tbin_start, tbin_stop = sim.getTimeSlices(GTI=(delay, tmax), return_bins=True) 

    # -------------------------------------------------------- simulate ---!!!
    print(f'Simulate template seed={sim.seed}')
    for j in range(tbin_stop-tbin_start-1):
        sim.t = [tgrid[j]+onset, tgrid[j + 1]+onset]
        if args.print:
            print(f'GTI (bin) = {sim.t} s')
        sim.model = join(datapath, f'extracted_data/{runid}/{runid}_tbin{tbin_start+j:02d}.xml')
        event = join(sim_output_path, f'{name}_tbin{tbin_start+j:02d}.fits')
        event_bins.append(event)
        sim.output = event
        try:
            sim.run_simulation()
        except Exception as e:
            rmtree(sim_output_path)
            return TrialOutput(trial_id, runid, time()-start_t, True, f"Simulation failed: {e}")

    # -------------------------------------------- shift time --- !!!
    if onset != 0:
        if delay != 0:
            return TrialOutput(trial_id, runid, time()-start_t, True, f'Bad configuration. Either "onset" or "delay" must be equal to 0.')
                
        # ------------------------------------ add background --- !!!
        print('Simulate bkg to append before the burst')
        bkg_name = f'bkg_{count:08d}_onset_{onset}_delay_{delay}_offset_{offset}'
        bkg = os.path.join(sim_output_path, bkg_name)
        event_bins.insert(0, bkg)
        sim.t = [0, onset]
        if args.print:
            print(f"GTI (bkg) = {sim.t} s")
        sim.model = bkg_model
        sim.output = bkg
        sim.run_simulation()

    # ---------------------------------------- gather bins ---!!!
    if args.merge:
        print('Merge in photon-list')
        phlist = join(sim_output_path, f'{name}.fits')
        sim.input = event_bins
        sim.output = phlist
        sim.appendEventsSinglePhList(GTI=[delay, delay+cfg.get('tobs')])
        if args.print:
            h = fits.open(phlist)
            print('Check GTI and EVENTS time range:')
            print('************')
            print(h[2].data)
            print(h[1].data.field('TIME').min(), h[1].data.field('TIME').max())
            print('************')
            h.close()
    else:
        # observation list ---!
        obslist = join(sim_output_path, f'{name}.xml')
        if os.path.isfile(obslist):
            os.remove(obslist)
        make_obslist(obslist=obslist, items=event_bins, names=name)

    sim.input = phlist
    sim.sortObsEvents()
    del sim

    # move final file and remove the temporary directory ---!
    if args.remove and args.merge:
        move(phlist, sim_output_path_base.joinpath(f'{name}.fits')) 
        rmtree(sim_output_path)

    # time ---!   
    elapsed_t = time()-start_t

    return TrialOutput(trial_id, runid, elapsed_t)
    
def yield_batch(batch_size, runids, trials, cfg, args):
    trials_batch = []
    for runid in runids:
        for trial in range(trials):
            trials_batch.append((runid, trial, cfg, args))
            if batch_size == len(trials_batch):
                yield trials_batch
                trials_batch = []
    if len(trials_batch) > 0:
        yield trials_batch

def stats(output_dir, trials_outputs):
    trials_elapsed_times_mean = np.array([trial_output.elapsed_time for trial_output in trials_outputs if not trial_output.errors]).mean()
    trials_elapsed_times_std = np.array([trial_output.elapsed_time for trial_output in trials_outputs if not trial_output.errors]).std()
    trials_with_errors = [trial_output for trial_output in trials_outputs if trial_output.errors]
    print(f"Trials elapsed time (mean): {trials_elapsed_times_mean} +- {trials_elapsed_times_std}\nNumber of trials with errors: {len(trials_with_errors)}", flush=True)
    with open(output_dir.joinpath("trials_with_errors.txt"), "a") as f:
        for trial_output in trials_with_errors:
            f.write(f"\n{trial_output}")

def main():
    
    parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
    parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
    parser.add_argument('--merge', type=str2bool, default=True, help='Merge in single phlist (true) or use observation library (false)')
    parser.add_argument('--remove', type=str2bool, default=True, help='Keep only outputs')
    parser.add_argument('--output-dir', type=str, default=None, help='Output directory')
    parser.add_argument('--print', type=str2bool, default=False, help='Print out results')
    parser.add_argument('-mpt', '--mp-threads', type=int, default=4, help='The size of the threads pool') 
    args = parser.parse_args()

    if args.remove and not args.merge:
        raise ValueError('Keyword "remove" cannot be True if keyword "merge" is False.')

    cfg = Config(args.cfgfile)
    runids = config_runid(cfg)
    trials = cfg.get('trials')

    args.output_dir = Path(args.output_dir)
    if not args.output_dir.exists():
        args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"root output dir: {args.output_dir}")
    print(f"total runids: {len(runids)}")
    print(f"total trials: {trials}")
    print(f"total simulations: {len(runids)*trials}")
    print(f"Threads: {args.mp_threads}")
    
    if cfg.get('simtype') == "bkg":
        simulate_trial = simulate_trial_bkg
    else:
        create_output_dirs(args.output_dir, cfg, runids)
        simulate_trial = simulate_trial_grb

    pool = Pool(args.mp_threads)

    count = 0
    for trial_batch in yield_batch(args.mp_threads, runids, trials, cfg, args):
        trials_outputs = pool.map(simulate_trial, trial_batch)
        count += len(trials_outputs)
        print(f"\n[{now()}] Batch processed, computing statistics..", flush=True)
        stats(args.output_dir, trials_outputs)
        print(f"Processed a total of {count} trials", flush=True)

    """Profiling run
    #import cProfile, pstats
    #profiler = cProfile.Profile()
    #profiler.enable()
    #trial_output = simulate_trial(('run0406_ID000126', 1, cfg, args))
    #profiler.disable()
    #stats = pstats.Stats(profiler).sort_stats('tottime')
    #stats.print_stats()"""
    
    print(f'\n...[{now()}] done.\n')


if __name__=='__main__':
    main()
