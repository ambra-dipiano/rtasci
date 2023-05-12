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
from os.path import isdir, join, isfile
from rtasci.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from rtasci.lib.RTAManageXml import ManageXml
from rtasci.lib.RTAUtils import *
from rtasci.cfg.Config import Config
from rtasci.aph.utils import *
from rtasci.lib.RTAUtilsGW import *

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
            logname = f"{datapath}/outputs/{runid}/{cfg.get('tool')}{cfg.get('type')}_offset{offset}_seed{start_count+1:06d}-{start_count+1+trials:06d}.txt"
            if isfile(logname):
                os.remove(logname)    # ------------------------------------------------------ loop trials ---!!!
            for i in range(trials):
                count = start_count + i + 1
                #print(f'seed = {count:06d}')
                name = f'ebl{count:06d}'
                if args.merge.lower() == 'true':
                    phlist = join(grbpath, name+'.fits')
                    model = phlist.replace('.fits', '_grb.xml').replace('/obs/', '/rta_products/')
                else:
                    phlist = join(grbpath, f'{name}.xml')
                    model = phlist.replace('.xml', '_grb.xml').replace('/obs/', '/rta_products/')
                os.system(f"cp {join(expandvars(cfg.get('model')), f'{runid}.xml')} {model}")
                fit = model.replace('grb.xml', 'fit.xml')
                if args.print.lower() == 'true':
                    print(f'Input observation: {phlist}')
                if not isfile(phlist):
                    print(f'Missing observation {phlist}. \nSkip runid {runid}.')
                    break

                # -------------------------------------------- loop exposure times ---!!!

                for exp in cfg.get('exposure'):   
                    if cfg.get('cumulative'):
                        times = increase_exposure(start=exp, stop=cfg.get('tobs'), function='linear')
                    if cfg.get('lightcurve'):
                        times = lightcurve_base_binning(start=cfg.get('delay'), stop=cfg.get('tobs'), exposure=exp)
                    else:
                        times = exp
                    if args.print.lower() == 'true':
                        print(f"Time selections = {times} s")

                    # ---------------------------------------------------------- loop binning ---!!!
                    for t in times:
                        if t == len(times) and cfg.get('lightcurve'):
                            break
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
                        if texp != exp:
                            raise ValueError("exp != texp")
                        if args.print.lower() == 'true':
                            print(f"Exposure = {texp} s")
                        grb.input = phlist
                        grb.output = selphlist
                        if args.merge.lower() == 'true':
                            grb.run_selection()
                        else:
                            prefix = join(grbpath, f'texp{exp}s_')
                            grb.run_selection(prefix=prefix)
                        # aperture photometry ---!
                        if '.fits' in selphlist:
                            results = photometrics_counts(selphlist, pointing=pointing, true_coords=target, events_type='events_filename')
                        elif '.xml' in selphlist:
                            results = photometrics_counts(selphlist, pointing=pointing, true_coords=target, events_type='events_list')
                        sigma = li_ma(results['on'], results['off'], results['alpha'])
                        if args.print.lower() == 'true':
                            print('Photometry counts:', results)
                            print('Li&Ma significance:', sigma)

                        # modify model
                        detection = ManageXml(model)
                        detection.modXml(overwrite=True)
                        detection.setTsTrue() 
                        detection.parametersFreeFixed(src_free=['Prefactor'])
                        detection.setModelParameters(parameters=['RA', 'DEC', 'Index'], values=[target[0], target[1], cfg.get('index')])
                        detection.closeXml()
                        # fit ---!
                        grb.input = selphlist
                        grb.model = model
                        grb.output = fit
                        grb.run_maxlikelihood()
                        # stats ---!
                        xml = ManageXml(fit)
                        try:
                            coords = xml.getRaDec()
                            ra = coords[0][0]
                            dec = coords[1][0]
                            ts = xml.getTs()[0]
                            sqrt_ts = np.sqrt(ts)
                            if args.print.lower() == 'true':
                                print('TS significance:', sqrt_ts)
                        except IndexError:
                            sqrt_ts = np.nan
                            print('Candidate not found.')
                        if sqrt_ts >= 0:
                            # flux ---!
                            spectra = xml.getSpectral()
                            index, pref, pivot = spectra[0][0], spectra[1][0], spectra[2][0]
                            err = xml.getPrefError()[0]
                            flux = phflux_powerlaw(index, pref, pivot, grb.e, unit='TeV')
                            flux_err = phflux_powerlaw(index, err, pivot, grb.e, unit='TeV')
                        else:
                            ra, dec, ts, sqrt_ts, flux, flux_err = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

                        if sqrt_ts < 5 and cfg.get('lightcurve') and cfg.get('cumulative'):
                            if exp < max(cfg.get('exposure')):
                                if args.print.lower() == 'true':
                                    print(f'No significance detection, increase exposure.')
                                break
                            else:
                                sys.exit(f"No significance detection with maximum exposure {texp} s.")
                                
                        if sigma < 5 or grb.t[1] > (cfg.get('tobs')+cfg.get('delay')):
                            break

                        row = f"{runid} {count} {grb.t[0]} {grb.t[1]} {exp} {sqrt_ts} {flux} {flux_err} {ra} {dec} {results['on']} {results['off']} {results['alpha']} {results['excess']} {sigma} {offset} {cfg.get('delay')} {cfg.get('scalefluxfactor')} {caldb} {irf} ctools3d_unbinned\n"
                        if args.print.lower() == 'true':
                            print(f"Results: {row}")
                        if not isfile(logname):
                            hdr = 'runid seed start stop texp sqrt_ts flux flux_err ra dec oncounts offcounts alpha excess sigma offset delay scaleflux caldb irf pipe\n'
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