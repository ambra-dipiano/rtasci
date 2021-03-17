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
from RTAscience.lib.RTACtoolsSimulation import RTACtoolsSimulation
from RTAscience.lib.RTAManageXml import ManageXml
from RTAscience.lib.RTAUtils import get_pointing, get_mergermap, get_alert_pointing_gw
from RTAscience.cfg.Config import Config
from os.path import isdir, isfile, join, expandvars

parser = argparse.ArgumentParser(description='Simulate empty fields.')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
args = parser.parse_args()
cfg = Config(args.cfgfile)

# general ---!
if cfg.get('simtype').lower() != 'bkg':
    raise ValueError('This script only allows bakground simulations')
trials = cfg.get('trials')  # trials
tobs = cfg.get('tobs')  # total obs time (s)

# paths ---!
datapath = cfg.get('data')
if not isdir(datapath):  # main data folder
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):  # obs parent folder
    os.mkdir(join(datapath, 'obs'))
bkgpath = join(datapath, 'obs', 'backgrounds')
if not isdir(bkgpath):
    os.mkdir(bkgpath)
# background model ---!
bkg_model = expandvars(cfg.get('bkg'))  # XML background model

# ------------------------------------------------------- loop runid --- !!!
pointing = [0., 0.]
if pointing[1] < 0:
    pointing[0] += 0.0
    pointing[1] += -cfg.get('offset')
else:
    pointing[0] += 0.0
    pointing[1] += cfg.get('offset')

for i in range(trials):
    count = cfg.get('start_count') + i + 1
    name = f'bkg{count:06d}'
    # setup ---!
    sim = RTACtoolsSimulation()
    sim.configure(cfg)
    sim.seed = count
    sim.pointing = pointing

    print('Simulate empty fields')
    sim.seed = count
    sim.t = [0, tobs]
    bkg = os.path.join(bkgpath, f'{name}.fits')
    sim.model = bkg_model
    sim.output = bkg
    sim.run_simulation()

print('.. done')