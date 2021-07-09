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
import numpy as np
from RTAscience.lib.RTACtoolsSimulation import RTACtoolsSimulation, make_obslist
from RTAscience.lib.RTAUtils import wobble_pointing
from RTAscience.cfg.Config import Config
from os.path import join, isfile, isdir, expandvars

parser = argparse.ArgumentParser(description='This scripts allows to simulate wobble observation')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
args = parser.parse_args()

cfg = Config(args.cfgfile)

# source ---!
if cfg.get('simtype').lower() != 'wobble':
    raise ValueError('Invalid simtype.')
runid = cfg.get('runid').lower()
if runid != 'crab':
    raise ValueError('Currenlty only "crab" runid are contemplated for wobble simulation.')

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
model = os.path.join(datapath, f'models/{runid}.xml')  # grb XML template model
obsfile = os.path.join(obspath, 'wobble.xml')

# ---------------------------------------------------- trials --- !!!
target = (83.63, 22.01)
runs = list()
for count in range(nruns):
    name = f'crab{start_count+count:06d}'
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
    run = os.path.join(obspath, f'{name}_run{count:02d}.fits')
    runs.append(run)
    sim.output = run
    sim.run_simulation()

make_obslist(obsfile, runs, runid)