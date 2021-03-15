# IMPORTS ---!
from pkg_blindsearch import *
import numpy as np
import os
import sys
import time
from os.path import isfile, isdir, join
import astropy.units as u
from astropy.coordinates import SkyCoord 
import argparse

# --------------------------------- SETUP --------------------------------- !!!

parser = argparse.ArgumentParser()
parser.add_argument('--trials', type=int, default=1000000, help='total trials')
parser.add_argument('--count', type=int, default=50000, help='trials per node')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--caldb', type=str, default='degr3b-v2', help='calibration database')
parser.add_argument('--irf', type=str, default='South_z40_0.5h', help='insturment response function')
args = parser.parse_args()

# initialize global count ---!
trials = args.trials  
start_count = args.count  
count = start_count  
offset = args.off
setdof = False 

# cpus ---!
nthreads = 1
os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
os.environ['MKL_NUM_THREADS'] = str(nthreads)

# ctools/cscripts parameters ---!
caldb = args.caldb
irf = args.irf

# fix all
texp = (10, 100)  # exposure times (s)
tint = len(texp)
tmin = 0  # slewing time (s)
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 150.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 150.0  # selection maximum energy (TeV)
roi = 5  # region of interest (deg)

ts_threshold = 25  # TS threshold for reliable detection

# conditions control ---!
checks = True  # prints checks info ---!
irf_degrade = True  # use degraded irf ---!
skip_exist = False  # skip trial if already existing in data file ---!
debug = False  # prints logfiles on terminal ---!
if_log = True  # saves logfiles ---!

# files ---!
model_pl = 'grb.xml'
model_bkg = 'CTAIrfBackground.xml'
tcsv = 'time_slices.csv'
merge_map = 'run0406_MergerID000126_skymap.fits'
ebl_template = 'run0406_ID000126_ebl.fits'
cfg = xmlConfig('/wilks.xml')
p = ConfigureXml(cfg)

# pointing with off-axis equal to max prob GW ---!
#pointing = (31.582, -53.211)  # pointing direction RA/DEC (deg)
true_coord, pointing, offmax = getPointing(fits_file=p.getWorkingDir()+ebl_template, merge_map=p.getWorkingDir()+merge_map)
if offset != 'gw':
    offset = float(offset)
    pointing = (true_coord[0], true_coord[1] - offset)
print('coords true:', true_coord, 'point', pointing, 'off', offmax) if checks else None


# recap and dof ---!
dof, m2, m1 = getDof()
print('!!! *** !!! dof = ', m2, ' - ', m1, ' = ', dof)
print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! caldb:', caldb)
print('!!! *** !!! irf:', irf)
print('!!! *** !!! sim energy range: [', elow, ', ', ehigh, '] (TeV)')
print('!!! *** !!! selection energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi: ', roi, ' (deg)')
print('!!! *** !!! pointing:', pointing, ' (deg)')

logname = f"{p.getCsvDir()}/{caldb}-{irf}_bkg{start_count+1:06d}-{start_count+trials:06d}_offset{offset}.txt"

if isfile(logname):
    os.remove(logname)

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup model dof ---!i
if setdof:
    model = p.getWorkingDir()+model_pl
    xml = ManageXml(model, '/wilks.xml')
    xml.prmsFreeFix()
    xml.closeXml()
    del xml
# setup trials obj ---!
bkg = Analysis('/wilks.xml')
bkg.nthreads = nthreads
bkg.pointing = pointing
bkg.roi = roi
bkg.e = [elow, ehigh]
bkg.tmax = max(tmax)
bkg.debug = debug
bkg.if_log = if_log
# degrade IRF if required ---!
bkg.caldb = caldb
bkg.irf = irf
print('!!! check ---- caldb:', bkg.caldb) if checks is True else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
    count += 1
    bkg.seed = count
    print('\n\n!!! ************ STARTING TRIAL %d ************ !!!' %count) if checks is True else None
    print('!!! check ---- seed=', bkg.seed) if checks is True else None

    # --------------------------------- SIMULATION --------------------------------- !!!

    # attach ID to fileroot ---!
    f = 'bkg%06d' % (count)
    # simulate ---!
    model = p.getWorkingDir() + model_bkg
    bkg.model = model
    event = p.getSimDir() + f + ".fits"
    if os.path.isfile(event):
        os.remove(event)
    bkg.output = event
    bkg.eventSim()
    print('!!! check ---- simulation=', event) if checks is True else None

    # --------------------------------- 2° LOOP :: texp --------------------------------- !!!

    # --------------------------------- SELECTION --------------------------------- !!!

    bkg.e = [emin, emax]
    for i in range(tint):
        print('\n\n!!! ************ STARTING TEXP %d ************ !!!\n\n' % texp[i]) if checks is True else None
        bkg.t = [tmin, tmax[i]]
        event_selected = event.replace(p.getSimDir(), p.getSelectDir()).replace('bkg%06d' %count, 'texp%ds_bkg%06d' %(texp[i], count))
        prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
        bkg.input = event
        bkg.output = event_selected
        bkg.eventSelect(prefix=prefix)
        print('!!! check ---- selection: ', event_selected) if checks is True else None

        # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

        model = p.getWorkingDir() + model_pl
        likeXml = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.fits', '_like.xml')
        if os.path.isfile(likeXml):
            os.remove(likeXml)
        bkg.input = event_selected
        bkg.model = model
        bkg.output = likeXml
        bkg.maxLikelihood()
        xml = ManageXml(likeXml, '/wilks.xml')
        print('!!! check ---- max likelihood: ', likeXml) if checks is True else None

        # --------------------------------- BEST FIT TSV --------------------------------- !!!

        ts_list, ts = ([] for j in range(2))
        ts_list.append(xml.loadTs())

        # only first elem ---!
        ts.append(ts_list[0][0])

        # --------------------------------- CLOSE LIKE XML ----------------- !!!

        xml.closeXml()

        # --------------------------------- RESULTS TABLE (csv) ----------------------- !!!

        if checks:
            print('\n!!! ---------- check trial:', count)
            print('!!! ----- check texp:', texp[i])
            print('!!! *** check ts:', ts[0])
        
        runid = 'run0406_ID000126'
        ts = float(ts[0])
        row = f"{runid} {count} {texp[i]} {ts} {offset} {caldb} {irf}\n"
        print(logname)
        print('!!! check row: seed %d --- texp' %i, texp[i], 's =====\n', row) if checks is True else None
        if not isfile(logname):
            hdr = 'runid seed texp ts offset caldb irf\n'
            log = open(logname, 'w+')
            log.write(hdr)
            log.write(row)
            log.close()
        else:
            log = open(logname, 'a')
            log.write(row)
            log.close()

    # --------------------------------- CLEAR SPACE --------------------------------- !!!

    print('!!! check ---- ', count, ') trial done...') if checks is True else None
    if int(count) != 20001:
        os.system('rm ' + p.getSimDir() + '*bkg%06d*' % count)
        os.system('rm ' + p.getSelectDir() + '*bkg%06d*' % count)
        os.system('rm ' + p.getDetDir() + '*bkg%06d*' % count)

print('\n\n!!! ================== END ================== !!!\n\n')
