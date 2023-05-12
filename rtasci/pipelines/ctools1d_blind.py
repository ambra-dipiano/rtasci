# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import time
tstamp = time.time()
import os
import sys
import argparse
import numpy as np
from os.path import isdir, join, isfile, expandvars
from rtasci.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from rtasci.lib.RTAManageXml import ManageXml
from rtasci.lib.RTAUtils import *
from rtasci.cfg.Config import Config
from rtasci.aph.utils import *
from rtasci.lib.RTAUtilsGW import *
from astropy.coordinates import SkyCoord

runtime = time.time() - tstamp

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

# CALDB ---!
if type(cfg.get('caldb')) == str:
    caldbs = [cfg.get('caldb')]
else:
    caldbs = cfg.get('caldb')
caldbs = sorted(caldbs)

# IRF ---!
if type(cfg.get('irf')) == str:
    irfs = [cfg.get('irf')]
else:
    irfs = cfg.get('irf')
irfs = sorted(irfs)

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
    if not isdir(f"{datapath}/outputs/{runid}"):
        os.mkdir(f"{datapath}/outputs/{runid}")
    if not isdir(f"{datapath}/rta_products/{runid}"):
        os.mkdir(f"{datapath}/rta_products/{runid}")
    png = f"{datapath}/skymaps/{runid}"
    if not isdir(png):
        os.mkdir(png)
    # grb path ---!
    grbpath = join(datapath, 'obs', runid)  
    if not isdir(grbpath):
        raise FileExistsError(f"Directory {runid} not found in {datapath}/obs")
    rtapath = f'{datapath}/rta_products/{runid}'
    # true coords ---!

    target = get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits")
    if args.print.lower() == 'true':
        print(f'Target true = {target} deg')
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
    if args.print.lower() == 'true':
        print(f'Pointing = {pointing} deg')

    # ------------------------------------------------------ loop caldb ---!!!
    for caldb in caldbs:
        if args.print.lower() == 'true':
            print(f'Calibration database: {caldb}')       
        # ------------------------------------------------------ loop irf ---!!!
        for irf in irfs:
            if args.print.lower() == 'true':
                print(f'Instrument response function: {irf}')  
            erange = check_energy_thresholds(erange=[cfg.get('emin'), cfg.get('emax')], irf=irf)
            # outputs
            logname = f"{datapath}/outputs/{runid}/{cfg.get('tool')}{cfg.get('type')}-{caldb}-{irf}-seed{start_count+1:06d}-{start_count+trials:06d}.txt"
            if isfile(logname):
                os.remove(logname)
            # ------------------------------------------------------ loop trials ---!!!
            for i in range(trials):
                count = start_count + i + 1
                if args.print.lower() == 'true':
                    print(f'Seed = {count:06d}')
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

                # --------------------------------------------------- loop exposure times ---!!!
                for exp in cfg.get('exposure'):   
                    if cfg.get('cumulative'):
                        times = increase_exposure(start=exp, stop=cfg.get('tobs'), function='linear')
                    elif cfg.get('lightcurve'):
                        print('here')
                        times = lightcurve_base_binning(start=cfg.get('delay'), stop=cfg.get('tobs'), exposure=exp)
                    else:
                        times = [exp]
                    if args.print.lower() == 'true':
                        print(f"Time selections = {times} s")

                    # ---------------------------------------------------------- loop binning ---!!!
                    for t in times:
                        tcpu = time.time()
                        if t == len(times) and cfg.get('lightcurve'):
                            break
                        # selection ---!
                        selphlist = phlist.replace(f'{name}', f'texp{exp}s_{name}')
                        grb = RTACtoolsAnalysis()
                        grb.caldb = caldb
                        grb.irf = irf
                        grb.roi = cfg.get('roi')
                        grb.e = erange
                        if cfg.get('lightcurve'):
                            grb.t = [t, t + exp]
                        else: 
                            grb.t = [cfg.get('delay'), cfg.get('delay')+t]
                        if args.print.lower() == 'true':
                            print(f"Selection t = {grb.t} s")
                        texp = grb.t[1] - grb.t[0]
                        if args.print.lower() == 'true':
                            print(f"Exposure = {texp} s")
                        grb.input = phlist
                        grb.output = selphlist
                        if args.merge.lower() == 'true':
                            grb.run_selection()
                        else:
                            prefix = join(grbpath, f'texp{exp}s_')
                            grb.run_selection(prefix=prefix)

                         # on/off ---!
                        if '.fits' in selphlist:
                            events_type = 'events_filename'
                        else:
                            events_type = 'events_list'
                            filenames = ManageXml(selphlist)
                            run_list = filenames.getRunList()
                            filenames.closeXml()
                            del filenames
                            selphlist = Photometrics.load_data_from_fits_file(run_list[0])
                            for file in run_list[1:]:
                                selphlist = np.append(selphlist, Photometrics.load_data_from_fits_file(file))

                        # on/off ---!
                        if '.fits' in selphlist:
                            onoff = selphlist.replace('.fits', '_cspha.xml').replace('/obs/', '/rta_products/')
                            events_type = 'events_filename'
                        else:
                            onoff = selphlist.replace('.xml', '_cspha.xml').replace('/obs/', '/rta_products/')
                            events_type = 'events_list'
                            filenames = ManageXml(selphlist)
                            run_list = filenames.getRunList()
                            filenames.closeXml()
                            del filenames
                            selphlist = Photometrics.load_data_from_fits_file(run_list[0])
                            for file in run_list[1:]:
                                selphlist = np.append(selphlist, Photometrics.load_data_from_fits_file(file))
                        
                        # skymap ---!
                        grb.input = selphlist
                        grb.output = sky
                        grb.run_skymap(wbin=cfg.get('skypix'), roi_factor=cfg.get('skyroifrac'))
                        # blind-search ---!
                        grb.sigma = cfg.get('sgmthresh')
                        grb.corr_rad = cfg.get('smooth')
                        grb.max_src = cfg.get('maxsrc')
                        grb.input = sky
                        grb.output = candidates
                        grb.run_blindsearch()
                        # position ---!
                        xml = ManageXml(candidates)
                        try:
                            target = xml.getRaDec()
                            ra = target[0][0]
                            dec = target[1][0]
                            offset = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs').separation(SkyCoord(ra=pointing[0], dec=pointing[1], unit='deg', frame='icrs')).deg
                            target = (ra, dec)
                            if args.print.lower() == 'true':
                                print(f'Target = [{ra}, {dec}]')
                            # on/off ---!
                            if '.fits' in selphlist:
                                onoff = selphlist.replace('.fits', '_cspha.xml').replace('/obs/', '/rta_products/')
                                events_type = 'events_filename'
                            else:
                                onoff = selphlist.replace('.xml', '_cspha.xml').replace('/obs/', '/rta_products/')
                                events_type = 'events_list'
                                filenames = ManageXml(selphlist)
                                run_list = filenames.getRunList()
                                filenames.closeXml()
                                del filenames
                                selphlist = Photometrics.load_data_from_fits_file(run_list[0])
                                for file in run_list[1:]:
                                    selphlist = np.append(selphlist, Photometrics.load_data_from_fits_file(file))

                            # aperture photometry ---!
                            phm = Photometrics({events_type: selphlist})
                            opts = phm_options(erange=grb.e, texp=texp, time_int=grb.t, target=target, pointing=pointing, index=-2.1, save_off_reg=f"{expandvars(cfg.get('data'))}/rta_products/{runid}/texp{exp}s_{name}_off_regions.reg", irf_file=join(expandvars('$CTOOLS'), f"share/caldb/data/cta/{caldb}/bcf/{irf}/irf_file.fits"))
                            off_regions = find_off_regions(phm, opts['background_method'], target, pointing, opts['region_radius'], verbose=opts['verbose'], save=opts['save_off_regions'])
                            on, off, alpha, excess, sigma, err_note = counting(phm, target, opts['region_radius'], off_regions, e_min=opts['energy_min'], e_max=opts['energy_max'], t_min=opts['begin_time'], t_max=opts['end_time'], draconian=False)
                            if args.print.lower() == 'true':
                                print(f'Photometry on={on} off={off} ex={excess} a={alpha}')
                                print('Li&Ma significance:', sigma)

                        except KeyError:
                            if args.print.lower() == 'true':
                                print('No candidates found.')   
                            ra, dec, on, off, alpha, excess, sigma, offset = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan , np.nan

                        if sigma < 5 or offset > cfg.get('roi')-0.5 or grb.t[1] > (cfg.get('tobs')+cfg.get('delay')):
                            break

                        else:
                            # flux ---!
                            src = {'ra': target[0], 'dec': target[1], 'rad': opts['region_radius']}
                            conf = ObjectConfig(opts)
                            region_eff_resp = aeff_eval(conf, src, {'ra': pointing[0], 'dec': pointing[1]})
                            livetime = opts['end_time'] - opts['begin_time']
                            flux = excess / region_eff_resp / livetime
                            gamma = opts['power_law_index']
                            if args.print.lower() == 'true':
                                print(f'Flux={flux}')

                        if sigma < 5 and cfg.get('lightcurve'):
                            if exp < max(cfg.get('exposure')) or grb.t[0] > cfg.get('tobs'):
                                if args.print.lower() == 'true':
                                    print(f'No significant detection, increase exposure.')
                                    break
                                else:
                                    sys.exit(f"No significant detection with max. exposure {texp} s.")

                        # save results ---!
                        k0, e0, flux_err, sqrt_ts = np.nan, np.nan, np.nan, np.nan
                        timing = runtime + time.time() - tcpu
                        row = f"{runid} {count} {grb.t[0]} {grb.t[1]} {grb.t[1]-grb.t[-0]} {sqrt_ts} {flux} {flux_err} {ra} {dec} {k0} {gamma} {e0} {on} {off} {alpha} {excess} {sigma} {offset} {cfg.get('delay')} {cfg.get('scalefluxfactor')} {caldb} {irf} ctools1d_blind {runtime}\n"
                        if args.print.lower() == 'true':
                            print(f"Results: {row}")
                        if not isfile(logname):
                            hdr = 'runid seed start stop texp sqrt_ts flux flux_err ra dec prefactor index scale on off alpha excess sigma offset delay scaleflux caldb irf pipe runtime\n'
                            log = open(logname, 'w+')
                            log.write(hdr)
                            log.write(row)
                            log.close()
                        else:
                            log = open(logname, 'a')
                            log.write(row)
                            log.close()

                        del grb, phm
                if args.remove.lower() == 'true':
                    # remove files ---!
                    os.system(f"rm {datapath}/obs/{runid}/texp*{name}*")
                    os.system(f"rm {datapath}/rta_products/{runid}/*{name}*")
print('...done.\n')


