

import os
import pandas as pd
from module_plot import *
from RTAscience.lib.RTACtoolsSimulation import RTACtoolsSimulation, RTACtoolsBase
from RTAscience.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from os.path import join

# files and path ---!
datapath = '/home/ambra/Desktop/CTA/projects/DATA/'
model = join(datapath, 'models/grb.xml')

cfg = xmlConfig('/config_irf.xml')
p = ConfigureXml(cfg=cfg)
path = p.getWorkingDir()
outpath = p.getRunDir()
pngpath = p.getPngDir()
# irf and caldb ---!
caldb_nom = 'prod3b-v2'
caldb_deg = caldb_nom.replace('prod', 'degr')
irf = 'South_z40_0.5h'
# setup ---!
nominal = True
degraded = True
compute = True
plot = False
sens_type = 'Integral'
print('Compute', sens_type, 'sensitivity')

caldb = []
if nominal:
  caldb.append(caldb_nom)
  print('Use nominal cladb')
if degraded:
  caldb.append(caldb_deg)
  print('Use degraded caldb')

e = [0.1, 0.44]
#texp = [42, 40, 70, 151.5, 438.5, 1654] # MAGIC 2
#texp = [38.2, 39.8, 70.28] # MAGIC 1
#texp = [1, 5, 10, 100] # CTA
texp = [400, 500, 600] # HESS
pointing = (33.057, -51.841) # pointing direction RA/DEC (deg) - centered on the source
nbins = 1 # energy bins for sensitivity computation
src_name = 'GRB'

# INITIALIZE ---!
if compute:
  for i in range(len(caldb)):
    for j in range(len(texp)):
      print('texp = ', texp[j], ' s')
      event = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_phlist.fits'
      results = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_maxlike.xml'
      output = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_HESSsens.csv'
      nObj = Analysis('/config_irf.xml')
      nObj.e = e
      nObj.t = [0, texp[j]]
      nObj.caldb = caldb[i]
      nObj.irf = irf
      nObj.model = model
      # NOMINAL SIM ---!
      nObj.output = event
      nObj.pointing = pointing  # (deg)
      nObj.eventSim()
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

    showSensitivity([10**energy_nom, 10**energy_deg], [sens_nom, sens_deg], savefig=savefig1[j], marker=['+', 'x'],
                    xlabel='energy (TeV)', ylabel='E$^2$F sensitivity (erg/cm$^2$/s)',
                    label=['full sens irf', 'degraded irf'], title=title, fontsize=12, show=False)

    showSensitivity([10**energy_nom, 10**energy_deg], [flux_nom, flux_deg], savefig=savefig2[j], marker=['+', 'x'],
                    xlabel='energy (TeV)', ylabel='ph flux (ph/cm$^2$/s)',
                    label=['full sens irf', 'degraded irf'], title=title, fontsize=12, show=False)