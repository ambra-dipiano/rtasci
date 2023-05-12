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
from os.path import isdir, isfile, join, expandvars
from rtasci.cfg.Config import Config
from rtasci.lib.RTAManageXml import ManageXml
from rtasci.lib.RTACtoolsSimulation import RTACtoolsSimulation
from rtasci.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from rtasci.lib.RTAUtils import get_pointing, get_mergermap
from rtasci.lib.RTAUtilsGW import get_alert_pointing_gw
from rtasci.lib.RTAVisualise import plotSkymap
from rtasci.aph.utils import *

def main(args):
    cfg = Config(args.cfgfile)
    # general ---!
    if cfg.get('simtype').lower() != 'bkg':
        raise ValueError('This script should be used only for empty fields simulation and analysis')
    trials = cfg.get('trials')  # trials
    tobs = cfg.get('tobs')  # total obs time (s)
    runid = cfg.get('runid')
    if type(runid) != str or runid == 'all':
        raise ValueError('Currently only 1 runid is allowed.')

    # paths ---!
    datapath = cfg.get('data')
    if not isdir(datapath):  # main data folder
        raise ValueError('Please specify a valid path')
    if not isdir(join(datapath, 'obs')):  # obs parent folder
        os.mkdir(join(datapath, 'obs'))
    bkgpath = join(datapath, 'obs', f'{runid}_backgrounds')
    if not isdir(bkgpath):
        os.mkdir(bkgpath)
    if not isdir(join(datapath, 'rta_products')):  
        os.mkdir(join(datapath, 'rta_products'))
    rtapath = join(datapath, 'rta_products', f'{runid}_backgrounds')
    if not isdir(rtapath):
        os.mkdir(rtapath)
    if not isdir(join(datapath, 'skymaps')): 
        os.mkdir(join(datapath, 'skymaps'))
    png = join(datapath, 'skymaps', f'{runid}_backgrounds')
    if not isdir(join(datapath, 'outputs')): 
        os.mkdir(join(datapath, 'outputs'))
    outpath = join(datapath, 'outputs', f'{runid}_backgrounds')
    if not isdir(outpath):
        os.mkdir(outpath)
    # background model ---!
    bkg_model = expandvars(cfg.get('bkg'))  # XML background model
    logname = join(outpath, f"{cfg.get('caldb')}-{cfg.get('irf')}_seed{cfg.get('start_count')+1:06d}-{cfg.get('start_count')+1+trials:06d}_offset{cfg.get('offset')}.txt")
    # true coords ---!
    true_coords = get_pointing(f"{os.path.expandvars(cfg.get('catalog'))}/{runid}.fits")

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
        
    for i in range(trials):
        # initialise ---!
        count = cfg.get('start_count') + i + 1
        name = f'bkg{count:06d}'
        print(f"Simulate empty fields for runid = {cfg.get('runid')} seed = {count}")

        # setup ---!
        sim = RTACtoolsSimulation()
        sim.seed = count
        sim.pointing = pointing
        sim.caldb = cfg.get('caldb')
        sim.irf = cfg.get('irf')
        sim.fov = cfg.get('roi')
        sim.e = [cfg.get('emin'), cfg.get('emax')]
        sim.t = [0, tobs]
        bkg = join(bkgpath, f'{name}.fits')
        sim.model = bkg_model
        sim.output = bkg
        sim.run_simulation()
        
        # analysis ---!
        for texp in cfg.get('exposure'):
            selphlist = bkg.replace(f'{name}', f'texp{texp}s_{name}')
            if args.print.lower() == 'true':
                print(f"Selection: {selphlist}")
            an = RTACtoolsAnalysis()
            an.pointing = pointing
            an.caldb = cfg.get('caldb')
            an.irf = cfg.get('irf')
            an.fov = cfg.get('roi')
            an.e = [cfg.get('emin'), cfg.get('emax')]
            an.t = [0, texp]
            an.input = bkg
            an.output = selphlist
            an.run_selection()
            # aperture photometry ---!
            if '.fits' in selphlist:
                results = photometrics_counts(selphlist, pointing=pointing, true_coords=true_coords, events_type='events_filename')
            elif '.xml' in selphlist:
                results = photometrics_counts(selphlist, pointing=pointing, true_coords=true_coords, events_type='events_list')
            sigma = li_ma(results['on'], results['off'], results['alpha'])
            if args.print.lower() == 'true':
                print('Photometry counts:', results)
                print('Li&Ma significance:', sigma)
            # skymap ---!
            sky = selphlist.replace(bkgpath, rtapath).replace('.fits', '_sky.fits')
            if args.print.lower() == 'true':
                print(f"Skymap: {sky}")
            an.input = selphlist
            an.output = sky
            an.max_src = cfg.get('maxsrc')
            an.sky_subtraction = 'NONE'
            an.run_skymap(wbin=cfg.get('skypix'), roi_factor=cfg.get('skyroifrac'))    
            # blind-search ---!
            candidates = sky.replace('_sky.fits', '_sources.xml')
            if args.print.lower() == 'true':
                print(f"Candidates: {candidates}")
            an.sigma = cfg.get('sgmthresh')
            an.corr_rad = cfg.get('smooth')
            an.max_src = cfg.get('maxsrc')
            an.input = sky
            an.output = candidates
            an.run_blindsearch()
            if cfg.get('plotsky'):
                plotSkymap(sky, reg=candidates.replace('.xml', '.reg'), suffix=f'{texp}s', png=png)
            # modify model
            detection = ManageXml(candidates)
            detection.modXml(overwrite=True)
            detection.setTsTrue() 
            detection.parametersFreeFixed(src_free=['Prefactor'])
            detection.closeXml()
            # fit ---!
            fit = candidates.replace('_sources.xml', '_fit.xml')
            if args.print.lower() == 'true':
                print(f"Fit: {fit}")
            an.input = selphlist
            an.model = candidates
            an.output = fit
            an.run_maxlikelihood()
            # stats ---!
            print(fit)
            xml = ManageXml(fit)
            try:
                coords = xml.getRaDec()
                ra = coords[0][0]
                dec = coords[1][0]
                ts = xml.getTs()[0]
            except IndexError:
                ts, ra, dec = np.nan, np.nan, np.nan
                print('Candidate not found.')

            row = f"{runid} {count} {texp} {ts} {ra} {dec} {results['on']} {results['off']} {results['alpha']} {results['excess']} {sigma} {cfg.get('offset')} {cfg.get('caldb')} {cfg.get('irf')}\n"
            if args.print.lower() == 'true':
                print(f"Results: {row}")
            if not isfile(logname):
                hdr = 'runid seed texp ts ra dec oncounts offcounts alpha excess sigma offset caldb irf\n'
                log = open(logname, 'w+')
                log.write(hdr)
                log.write(row)
                log.close()
            else:
                log = open(logname, 'a')
                log.write(row)
                log.close()

            del an
        del sim
        if args.remove.lower() == 'true':
            # remove files ---!
            os.system(f"rm {datapath}/obs/{runid}_backgrounds/*{name}*")
            os.system(f"rm {datapath}/rta_products/{runid}_backgrounds/*{name}*")



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Simulate empty fields.')
    parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
    parser.add_argument('--remove', type=str, default='true', help='Keep only outputs')
    parser.add_argument('--print', type=str, default='false', help='Print out results')
    args = parser.parse_args()

    print(args.cfgfile)

    main(args)
    print('...done.\n')

