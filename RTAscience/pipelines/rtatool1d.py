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
from os.path import isdir, join, isfile, expandvars
from RTAscience.lib.RTACtoolsAnalysis import RTACtoolsAnalysis, onoff_counts
from RTAscience.lib.RTAManageXml import ManageXml
from RTAscience.lib.RTAUtils import phflux_powerlaw, get_pointing, get_mergermap, get_alert_pointing_gw, phmOptions, increase_exposure
from RTAscience.cfg.Config import Config
from RTAscience.lib.RTAVisualise import plotSkymap
from RTAscience.aph.utils import *

parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
parser.add_argument('--merge', type=str, default='true', help='Merge in single phlist (true) or use observation library (false)')
parser.add_argument('--remove', type=str, default='true', help='Keep only outputs')
parser.add_argument('--print', type=str, default='false', help='Print out results')
args = parser.parse_args()

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
start_count = cfg.get('start_count')
trials = cfg.get('trials') 
if cfg.get('offset') == 'str':
    offset = cfg.get('offset').upper()
else:
    offset = cfg.get('offset')
# paths ---!
datapath = cfg.get('data')
if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):  # obs parent folder
    raise ValueError(f'Missing obs parent folder in {datapath}')
if not isdir(f"{datapath}/outputs"):
    os.mkdir(f"{datapath}/outputs")
if not isdir(f"{datapath}/rta_products"):
    os.mkdir(f"{datapath}/rta_products")
if not isdir(f"{datapath}/skymaps"):
    os.mkdir(f"{datapath}/skymaps")

# ------------------------------------------------------ loop runid --- !!!
for runid in runids:
    print(f"{'-'*50} #\nProcessing runid: {runid}")
    # outputs
    logname = f"{datapath}/outputs/{runid}/{cfg.get('caldb')}-{cfg.get('irf')}_seed{start_count+1:06d}-{start_count+1+trials:06d}_flux{cfg.get('scalefluxfactor')}_offset{offset}_delay{cfg.get('delay')}.txt"
    if not isdir(f"{datapath}/outputs/{runid}"):
        os.mkdir(f"{datapath}/outputs/{runid}")
    if not isdir(f"{datapath}/rta_products/{runid}"):
        os.mkdir(f"{datapath}/rta_products/{runid}")
    png = f"{datapath}/skymaps/{runid}"
    if not isdir(png):
        os.mkdir(png)
    if isfile(logname):
        os.remove(logname)
    # grb path ---!
    grbpath = join(datapath, 'obs', runid)  
    if not isdir(grbpath):
        raise FileExistsError(f"Directory {runid} not found in {datapath}/obs")
    rtapath = f'{datapath}/rta_products/{runid}'
    # true coords ---!

    target = get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits")
    # get alert pointing
    if type(cfg.get('offset')) == str and cfg.get('offset').lower() == 'gw':
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

    # ------------------------------------------------------ loop trials ---!!!
    for i in range(trials):
        count = start_count + i + 1
        #print(f'seed = {count:06d}')
        name = f'ebl{count:06d}'
        if args.merge.lower() == 'true':
            phlist = join(grbpath, name+'.fits')
            sky = phlist.replace('.fits', '_sky.fits').replace('/obs/', '/rta_products/')
        else:
            phlist = join(grbpath, f'{name}.xml')
            sky = phlist.replace('.xml', '_sky.fits').replace('/obs/', '/rta_products/')
        candidates = sky.replace('_sky.fits', '_sources.xml')
        fit = candidates.replace('sources', 'fit')
        if args.print.lower() == 'true':
            print(f'Input observation: {phlist}')
        if not isfile(phlist):
            print(f'Missing observation {phlist}. \nSkip runid {runid}.')
            break

        # ---------------------------------------------------------- loop exposure times ---!!!

        if cfg.get('cumulative'):
            n = int(cfg.get('tobs') / cfg.get('exposure')[0])
            times = [cfg.get('exposure')[0]*(i+1) for i in range(n)]
            if times[-1] < cfg.get('tobs'):
                times.append(cfg.get('tobs'))
        else:
            times = cfg.get('exposure')
        if args.print.lower() == 'true':
            print(f"Time selections = {times} s")
        # selection ---!
        for texp in times:
            if args.print.lower() == 'true':
                print(f"Exposure = {texp} s")
            selphlist = phlist.replace(f'{name}', f'texp{texp}s_{name}')
            grb = RTACtoolsAnalysis()
            grb.caldb = cfg.get('caldb')
            grb.irf = cfg.get('irf')
            grb.roi = cfg.get('roi')
            grb.e = [cfg.get('emin'), cfg.get('emax')]
            grb.t = [cfg.get('delay'), cfg.get('delay')+texp]
            if args.print.lower() == 'true':
                print(f"Selection t = {grb.t} s")
            grb.input = phlist
            grb.output = selphlist
            if args.merge.lower() == 'true':
                grb.run_selection()
            else:
                prefix = join(grbpath, f'texp{texp}s_')
                grb.run_selection(prefix=prefix)

            # on/off ---!
            if '.fits' in selphlist:
                onoff = selphlist.replace('.fits', '_cspha.xml').replace('/obs/', '/rta_products/')
                events_type = 'events_filename'
            else:
                onoff = selphlist.replace('.xml', '_cspha.xml').replace('/obs/', '/rta_products/')
                events_type = 'events_list'
                filenames = ManageXml(event_selected)
                run_list = filenames.getRunList()
                filenames.closeXml()
                del filenames
                selphlist = Photometrics.load_data_from_fits_file(run_list[0])
                for file in run_list[1:]:
                    selphlist = np.append(selphlist, Photometrics.load_data_from_fits_file(file))
            
            # aperture photometry ---!
            phm = Photometrics({events_type: selphlist})
            opts = phmOptions(cfg, texp=texp, target=target, pointing=pointing, runid=runid)
            off_regions = find_off_regions(phm, opts['background_method'], target, pointing, opts['region_radius'], verbose=opts['verbose'], save=opts['save_off_regions'])
            oncounts, offcounts, alpha, excess, sigma, err_note = counting(phm, target, opts['region_radius'], off_regions, e_min=opts['energy_min'], e_max=opts['energy_max'], t_min=opts['begin_time'], t_max=opts['end_time'], draconian=False)
            if args.print.lower() == 'true':
                print(f'Photometry on={oncounts} off={offcounts} ex={excess} a={alpha}')
                print('Li&Ma significance:', sigma)

            # flux ---!
            src = {'ra': target[0], 'dec': target[1], 'rad': opts['region_radius']}
            conf = ObjectConfig(opts)
            region_eff_resp = aeff_eval(conf, src, {'ra': pointing[0], 'dec': pointing[1]})
            livetime = opts['end_time'] - opts['begin_time']
            flux = excess / region_eff_resp / livetime
            k0, e0, flux_err, sqrt_ts = np.nan, np.nan, np.nan, np.nan
            ra = target[0]
            dec = target[1]
            gamma = opts['power_law_index']
            if args.print.lower() == 'true':
                print(f'Flux={flux}')

            # save results ---!
            row = f"{runid} {count} {texp} {sqrt_ts} {flux} {flux_err} {ra} {dec} {k0} {gamma} {e0} {oncounts} {offcounts} {alpha} {excess} {sigma} {offset} {cfg.get('delay')} {cfg.get('scalefluxfactor')} {cfg.get('caldb')} {cfg.get('irf')}\n"
            if args.print.lower() == 'true':
                print(f"Results: {row}")
            if not isfile(logname):
                hdr = 'runid seed texp sqrt_ts flux flux_err ra dec prefactor index scale oncounts offcounts alpha excess sigma offset delay scaleflux caldb irf\n'
                log = open(logname, 'w+')
                log.write(hdr)
                log.write(row)
                log.close()
            else:
                log = open(logname, 'a')
                log.write(row)
                log.close()

            del grb
        if args.remove.lower() == 'true':
            # remove files ---!
            os.system(f"rm {datapath}/obs/{runid}/texp*{name}*")
            os.system(f"rm {datapath}/rta_products/{runid}/*{name}*")
print('...done.\n')


