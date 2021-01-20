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
from lib.RTACtoolsSimulation import RTACtoolsSimulation
from lib.RTAUtils import wobble_pointing
from lib.RTACtoolsAnalysis import make_obslist

# source ---!
runid = 'crab'
# general ---!
trials = 4  # trials
count = 0  # starting count
nthreads = 2
# sim parameters ---!
caldb = 'prod3b'  # calibration database
irf = 'South_z20_average_LST_30m'  # istrument response function
nruns = 4  # total number of runs
trun = 1e2  # single run duration (s)
emin = 3e-2  # simulation minimum energy (TeV)
emax = 30.  # simulation maximum energy (TeV)
roi = 2.5  # region of interest radius (deg)
seed = 1  # seed of MC simulation
# paths ---!
pypath = str(os.path.dirname(os.path.abspath(__file__)))  
datapath = pypath.replace('cta-sag-sci', 'DATA')  # all data should be under this folder
obspath = os.path.join(datapath, 'obs', runid)  # folder that will host the phlist src+bkg phlists
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
    name = f'crab{count:06d}'
    wobble_index = count % 4
    pointing = wobble_pointing(target, nrun=count, clockwise=True, offset=0.5)
    print(f'Simulating run {count} at pointing: {pointing} deg')

    # setup ---!
    sim = RTACtoolsSimulation()
    sim.seed = seed
    sim.nthreads = nthreads
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