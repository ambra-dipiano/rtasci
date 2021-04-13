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

e = [0.03, 150.0]
#texp = [42, 40, 70, 151.5, 438.5, 1654] # MAGIC 2
#texp = [38.2, 39.8, 70.28] # MAGIC 1
texp = [1, 5, 10, 100] # CTA
#texp = [400, 500, 600] # HESS
pointing = (33.057, -51.841) # pointing direction RA/DEC (deg) - centered on the source
nbins = 0 # energy bins for sensitivity computation
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
      nObj.eventSens(bins=nbins)

# ------------------------------------- PLOT --------------------------------------- !!!

if plot:
  csv = [[], []]
  savefig1, savefig2, savefig3 = [], [], []
  list_sens_nom, list_flux_nom, list_sens_deg, list_flux_deg = [], [], [], []
  for i in range(len(caldb)):
    for j in range(len(texp)):
      csv[i].append(outpath + 'texp%ds_' %texp[j] + caldb[i] + '_crab_sens.csv')
      pngroot = caldb[i] + '_texp%ds' %texp[j]
      if sens_type.capitalize() != 'Integral':
        savefig1.append(pngpath + pngroot + '_sensDiff.png')
        savefig2.append(pngpath + pngroot + '_sensDiff_phflux.png')
      else:
        savefig1.append(pngpath + pngroot + '_sensInt.png')
        savefig2.append(pngpath + pngroot + '_sensInt_phflux.png')

  for j in range(len(texp)):
    title = caldb_nom + ': ' + irf.replace('_', '\_') + ' with texp=%ds' %texp[j]
    # nominal
    df_nom = pd.read_csv(csv[0][j])
    cols = list(df_nom.columns)
    energy_nom = np.array(df_nom[cols[0]])
    sens_nom = np.array(df_nom[cols[6]])
    flux_nom = np.array(df_nom[cols[4]])
    # degraded ---!
    df_deg = pd.read_csv(csv[1][j])
    energy_deg = np.array(df_deg[cols[0]])
    sens_deg = np.array(df_deg[cols[6]])
    flux_deg = np.array(df_deg[cols[4]])

