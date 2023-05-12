# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import sys
from rtasci.thesis.pkg_blindsearch import *
from os.path import join, expandvars, isdir

inst = str(sys.argv[1])

# files and path ---!
datapath = expandvars('$DATA')
model = join(datapath, 'models/grb.xml')

# files and path ---!
cfg = xmlConfig('/config.xml')
p = ConfigureXml(cfg=cfg)
outpath = join(datapath, 'irf_degradation/')
if not isdir(outpath):
    os.mkdir(outpath)
pngpath = join(outpath, 'png')
if not isdir(pngpath):
    os.mkdir(pngpath)
# irf and caldb ---!
caldb_nom = 'prod3b-v2'
caldb_deg = caldb_nom.replace('prod', 'degr')
irf = 'South_z40_0.5h'
p.setRunDir(irf + '/')
# setup ---!
nominal = True
degraded = True
compute = True
plot = True
sens_type = 'Integral'
print('Compute', sens_type, 'sensitivity')

caldb = []
if nominal:
  caldb.append(caldb_nom)
  print('Use nominal cladb')
if degraded:
  caldb.append(caldb_deg)
  print('Use degraded caldb')

if inst.lower() == 'magic':
  e = [0.3, 1.]
  #texp = [42, 40, 70, 151.5, 438.5, 1654] # MAGIC 2
  texp = [38.2, 39.8, 70.28] # MAGIC 1
elif inst.lower() == 'hess':
  e = [0.1, 0.44]
  texp = [300, 400, 500, 600] # HESS
elif inst.lower() == 'cta':
  e = [0.03, 150.0]
  #texp = [5, 10, 50, 100, 200, 300, 500] # CTA
  texp = [10, 100]

target = (33.057, -51.841) # centered on the source
offsets = [0]#, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
nbins = 10   # energy bins for sensitivity computation
src_name = 'GRB'

# INITIALIZE ---!
if compute:
  for i in range(len(caldb)):
    for offset in offsets:
      pointing = (target[0], target[1]-offset) 
      for j in range(len(texp)):
        event = outpath + 'texp%ds_' %texp[j] + caldb[i] + f'_{inst.upper()}sim.fits'
        results = outpath + 'texp%ds_' %texp[j] + caldb[i] + f'_{inst.upper()}fit.xml'
        output = outpath + 'texp%ds_' %texp[j] + caldb[i] + f'_{inst.upper()}sens.csv'
        nObj = Analysis('/config.xml')
        nObj.e = e
        nObj.t = [0, texp[j]]
        nObj.caldb = caldb[i]
        nObj.roi = 5
        nObj.irf = irf
        nObj.model = model
        # NOMINAL SIM ---!
        nObj.output = event
        nObj.pointing = pointing  # (deg)
        nObj.eventSim()
        print('caldb = ', caldb[i])
        print('texp = ', texp[j], ' s')
        print('off = ', offset)
        # NOMINAL MAX LIKELIHOOD ---!
        nObj.input = event
        nObj.output = results
        nObj.maxLikelihood()
        # NOMINAL SENS ---!
        nObj.sens_type = sens_type
        nObj.model = results
        nObj.output = output
        nObj.src_name = src_name
        nObj.eventSens(bins=nbins)

        if not isdir(f"{outpath}/off{offset}"):
          os.mkdir(f"{outpath}/off{offset}")
        os.system(f"mv {outpath}/*csv {outpath}/off{offset}/.")
