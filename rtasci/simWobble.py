# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import argparse
from rtasci.lib.RTACtoolsSimulation import RTACtoolsSimulation, make_obslist
from rtasci.lib.RTAUtils import wobble_pointing
from rtasci.cfg.Config import Config
from os.path import join, isdir

parser = argparse.ArgumentParser(description='This scripts allows to simulate wobble observation')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
args = parser.parse_args()

cfg = Config(args.cfgfile)

# source ---!
if cfg.get('simtype').lower() != 'wobble':
    raise ValueError('Invalid simtype.')
runid = cfg.get('runid').lower()
if runid not in ['crab', 'bkg']:
    raise ValueError('Currenlty only "crab" and "bkg" runid are contemplated for wobble simulation.')

# general ---!
trials = cfg.get('trials')
start_count = cfg.get('start_count') 
# sim parameters ---!
caldb = cfg.get('caldb')
irf = cfg.get('irf')
nruns = cfg.get('nruns')
trun = cfg.get('tobs')
emin = cfg.get('emin')
emax = cfg.get('emax')
roi = cfg.get('roi')
seed = start_count+1 

# paths ---!
datapath = cfg.get('data')
if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):  # obs parent folder
    os.mkdir(join(datapath, 'obs'))
obspath = os.path.join(datapath, 'obs', runid)  

# check folders and create missing ones ---!
if not os.path.isdir(obspath):
    if not os.path.isdir(obspath.replace(runid, '')):
        os.mkdir(obspath.replace(runid, ''))
    os.mkdir(obspath)
# files ---!
if runid == 'bkg':
    model = os.path.join(datapath, f'models/CTAIrfBackground.xml')
else:
    model = os.path.join(datapath, f'models/{runid}.xml')  # grb XML template model
obsfile = os.path.join(obspath, 'wobble.xml')

# ---------------------------------------------------- trials --- !!!
target = (83.6331, 22.0145)
runs = list()
for count in range(nruns):
    name = f'{cfg.get("runid")}{start_count:06d}'
    print(f'Name: {name}')
    wobble_index = count % 4
    pointing = wobble_pointing(target, nrun=count, clockwise=True, offset=cfg.get('offset'))
    print(f'Simulating run {count+1}/{nruns} at pointing: {pointing} deg')

    # setup ---!
    sim = RTACtoolsSimulation()
    sim.seed = seed
    sim.pointing = pointing
    sim.roi = roi
    sim.e = [emin, emax]
    sim.t = [trun*count, trun*(count+1)]
    sim.caldb = caldb
    sim.irf = irf
    sim.model = model
    run = os.path.join(obspath, f'{name}_run{count+1:02d}.fits')
    runs.append(run)
    sim.output = run
    sim.run_simulation()
    sim.input = run
    sim.sortObsEvents()

make_obslist(obsfile, runs, runid)