import os
import sys
import pandas as pd
from RTAscience.thesis.pkg_blindsearch import *
from os.path import join, expandvars, isdir, isfile

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
plot = False
sens_type = 'integral'

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
  texp = [5, 10, 100] # CTA

target = (33.057, -51.841) # centered on the source
offset = 1.633
pointing = (target[0], target[1]-offset) # with offset
nbins = 1   # energy bins for sensitivity computation
src_name = 'GRB'

# INITIALIZE ---!
if compute:
  for i in range(len(caldb)):
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
      print('texp = ', texp[j], ' s')
      print('caldb = ', caldb[i])
      # NOMINAL MAX LIKELIHOOD ---!
      nObj.input = event
      nObj.output = results
      nObj.maxLikelihood()
      # NOMINAL SENS ---!
      nObj.sens_type = sens_type
      nObj.model = results
      nObj.output = output
      nObj.src_name = src_name
      nObj.eventSens(bins=1)

os.system(f"mv {directory}/*csv {directory}/off{offset}")