# IMPORTS ---!
from pkg_blindsearch import *
import numpy as np
import csv
import os
import sys
import time
from os.path import isfile, isdir, join
import astropy.units as u
from astropy.coordinates import SkyCoord 
import argparse

# --------------------------------- SETUP --------------------------------- !!!

parser = argparse.ArgumentParser()
parser.add_argument('--trials', type=int, default=10000, help='total trials')
parser.add_argument('--count', type=int, default=500, help='trials per node')
parser.add_argument('--delay', type=float, default=50, help='delay')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--flux', type=float, default=1, help='flux scaling factor')
args = parser.parse_args()

# initialize global count ---!
trials = args.trials  
start_count = args.count  
count = start_count  
delay = args.delay 
offset = args.off
scaleflux = args.flux
extract_spec = False 

# cpus ---!
nthreads = 1
os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
os.environ['MKL_NUM_THREADS'] = str(nthreads)

# ctools/cscripts parameters ---!
caldb = 'degr3b-v2'
irf = 'South_z40_0.5h'

# fix all
sigma = 5  # detection acceptance (Gaussian)
texp = (10, 100)  # exposure times (s)
tint = len(texp)
tmin = [delay for i in range(tint)]  # start of bin to select (s)
tmax = []  # end of bin to select (s)
for i in range(tint):
    tmax.append(tmin[i] + texp[i])
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 150.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 150.0  # selection maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)
wbin = 0.02  # skymap bin width (deg)
corr_rad = 0.1  # Gaussian
confidence = (0.68, 0.95, 0.9973)  # confidence interval for asymmetrical errors (%)
max_src = 1  # max candidates
ts_threshold = 25  # TS threshold for reliable detection

# conditions control ---!
checks = True  # prints checks info ---!
if_ebl = True  # uses the EBL absorbed template ---!
if_cut = False  # adds a cut-off parameter to the source model ---!
ebl_fits = False  # generate the EBL absorbed template ---!
irf_degrade = True  # use degraded irf ---!
src_sort = True  # sorts scandidates from highest TS to lowest ---!
skip_exist = False  # skip trial if already existing in data file ---!
debug = False  # prints logfiles on terminal ---!
if_log = True  # saves logfiles ---!

# path configuration ---!
cfg = xmlConfig()
p = ConfigureXml(cfg)
# files ---!
fileroot = 'run0406_'
ebl_table = p.getRootDir() + 'ebl_tables/gilmore_tau_fiducial.csv'
merge_map = 'run0406_MergerID000126_skymap.fits'
nominal_template = 'run0406_ID000126.fits'
ebl_template = 'run0406_ID000126_ebl.fits'
model_pl = 'run0406_ID000126.xml'
tcsv = 'time_slices.csv'

# pointing with off-axis equal to max prob GW ---!
true_coord, pointing, offmax = getPointing(fits_file=p.getWorkingDir()+ebl_template, merge_map=p.getWorkingDir()+merge_map)
if offset != 'gw':
    offset = float(offset)
    pointing = (true_coord[0], true_coord[1] - offset)
print('coords true:', true_coord, 'point', pointing, 'off', offmax) if checks else None

# recap and dof ---!
dof, m2, m1 = getDof()
print('!!! *** !!! dof = ', m2, ' - ', m1, ' = ', dof)
print('!!! *** !!! EBL ABSORPTION:', if_ebl)
print('!!! *** !!! MODEL CUTOFF:', if_cut)
print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! nominal caldb:', caldb)
print('!!! *** !!! irf:', irf)
print('!!! *** !!! TS SORT:', src_sort)
print('!!! *** !!! FLUX REDUCED factor:', scaleflux)
print('!!! *** !!! sim energy range: [', elow, ', ', ehigh, '] (TeV)')
print('!!! *** !!! selection energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi: ', roi, ' (deg)')
print('!!! *** !!! target:', true_coord, '(deg)')
print('!!! *** !!! pointing:', pointing, ' (deg)')

logname = f"{p.getCsvDir()}/{caldb}-{irf}_seed{start_count+1:06d}-{start_count+trials:06d}_flux{scaleflux}_offset{offset}_delay{delay}.txt"

if isfile(logname):
    os.remove(logname)

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup trials obj ---!
grb = Analysis()
grb.nthreads = nthreads
grb.pointing = pointing
grb.roi = roi
grb.e = [elow, ehigh]
grb.tmax = max(tmax)
grb.model = p.getWorkingDir() + model_pl
grb.debug = debug
grb.if_log = if_log
# degrade IRF if required ---!
if irf_degrade:
    grb.caldb = caldb.replace('prod', 'degr')
else:
    grb.caldb = caldb
grb.irf = irf
# add EBL to template ---!
if ebl_fits:
    grb.template = p.getWorkingDir() + nominal_template # nominal ---!
    new_template = p.getWorkingDir() + ebl_template # absorbed ---!
    grb.table = ebl_table # fiducial ---!
    grb.zfetch = True
    grb.if_ebl = False
    grb.fitsEbl(new_template)
# assign template ---!
if if_ebl:
    template = p.getWorkingDir() + ebl_template
else :
    template = p.getWorkingDir() + nominal_template
grb.if_ebl = if_ebl
grb.template = template
# load template ---!
grb.extract_spec = extract_spec
if extract_spec:
    tbin_stop = grb.loadTemplate()
else:
    tbin_stop = grb.getTimeBinStop()

# --------------------------------- REDUCE TEMPLATE FLUX  --------------------------------- !!!

if scaleflux != 1:
    grb.factor = scaleflux
    if extract_spec:
        grb.makeFainter()
    print('!!! check ---- reduce flux by factor %s' %str(scaleflux)) if checks is True else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
    count += 1
    grb.seed = count
    print('!!! ************ STARTING TRIAL %d ************ !!!\n\n' %count) if checks is True else None
    print('!!! check ---- seed=', grb.seed) if checks is True else None
    # attach ID to fileroot ---!
    f = 'ebl%06d' % (count)
    print('!!! check ---- file=', f) if checks is True else None

    # --------------------------------- SIMULATION --------------------------------- !!!

    event_bins = []
    grb.table = p.getDataDir() + tcsv
    tgrid, start, stop = grb.getTimeSlicesNew(GTI=[delay, max(tmax)], return_bins=True)  
    # simulate ---!
    for i in range(len(tgrid)-1):
        grb.t = [tgrid[i], tgrid[i+1]]
        print(f'GTI={grb.t}')
        if if_ebl:
            grb.model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % (i+start)
            event = p.getSimDir() + f + "_ebl_tbin%02d.fits" % (i+start)
        else:
            grb.model = p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % (i+start)
            event = p.getSimDir() + f + "_tbin%02d.fits" % (i+start)
        if scaleflux != 1:
            grb.model = grb.model.replace('_tbin', '_flux%s_tbin' %str(scaleflux))
            event = event .replace('_tbin', '_flux%s_tbin' %str(scaleflux))
        event_bins.append(event)
        if os.path.isfile(event):
            os.remove(event)
        grb.output = event
        grb.eventSim()
    # observation list ---!
    event = event_bins
    event_list = p.getSimDir() + 'obs_%s.xml' % f
    if scaleflux != 1:
        event_list = event_list.replace('obs_', 'obs_flux%s_' %str(scaleflux))
    if os.path.isfile(event_list):
        os.remove(event_list)
    grb.input = event
    grb.output = event_list
    grb.obsList(obsname=f)

    # --------------------------------- 2° LOOP :: texp --------------------------------- !!!

    # --------------------------------- SELECTION --------------------------------- !!!

    grb.e = [emin, emax]
    for i in range(tint):
        print('!!! ************ STARTING TEXP %d ************ !!!\n\n' % texp[i]) if checks is True else None
        grb.t = [tmin[i], tmax[i]]
        event_selected = event_list.replace(p.getSimDir(), p.getSelectDir()).replace('obs_', 'texp%ds_' % texp[i])
        prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
        if os.path.isfile(event_selected):
            os.remove(event_selected)
        grb.input = event_list
        grb.output = event_selected
        grb.eventSelect(prefix=prefix)

        # --------------------------------- SKYMAP --------------------------------- !!!

        skymap = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.xml', '_skymap.fits')
        if os.path.isfile(skymap):
            os.remove(skymap)
        grb.input = event_selected
        grb.output = skymap
        grb.eventSkymap(wbin=wbin)

        # --------------------------------- DETECTION & MODELING --------------------------------- !!!

        grb.corr_rad = corr_rad
        grb.max_src = max_src
        detectionXml = skymap.replace('_skymap.fits', '_det%dsgm.xml' %sigma)
        if os.path.isfile(detectionXml):
            os.remove(detectionXml)
        grb.input = skymap
        grb.output = detectionXml
        grb.runDetection()
        degrb = ManageXml(detectionXml)
        degrb.sigma = sigma
        degrb.if_cut = if_cut
        degrb.modXml()
        degrb.prmsFreeFix()

        # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

        start_time = time.time()
        likeXml = detectionXml.replace('_det%dsgm' % grb.sigma, '_like%dsgm' % grb.sigma)
        if os.path.isfile(likeXml):
            os.remove(likeXml)
        grb.input = event_selected
        grb.model = detectionXml
        grb.output = likeXml
        grb.maxLikelihood()
        likeObj = ManageXml(likeXml)
        elapsed = time.time() - start_time
        print(elapsed, 's for %d' %texp[i])
        if src_sort:
            highest_ts_src = likeObj.sortSrcTs()[0]
        else:
            highest_ts_src = None
        print('!!! check ---- highest TS: ', highest_ts_src) if checks is True else None

        # --------------------------------- DETECTION RA & DEC --------------------------------- !!!

        pos, ra_det, dec_det = ([] for j in range(3))
        pos.append(degrb.loadRaDec(highest=highest_ts_src))
        ra_det.append(pos[0][0][0]) if len(pos[0][0]) > 0 else ra_det.append(np.nan)
        dec_det.append(pos[0][1][0]) if len(pos[0][0]) > 0 else dec_det.append(np.nan)
        Ndet = len(pos[0][0])

        # --------------------------------- CLOSE DET XML --------------------------------- !!!

        degrb.closeXml()

        # --------------------------------- BEST FIT TSV --------------------------------- !!!

        ts_list, ts = ([] for j in range(2))
        ts_list.append(likeObj.loadTs()) if Ndet > 0 else ts_list.append([np.nan])

        # only first elem ---!
        ts.append(ts_list[0][0])

        # --------------------------------- Nsrc FOR TSV THRESHOLD --------------------------------- !!!

        # count src with TS >= 9
        n = 0
        for j in range(len(ts_list[0])):
            if float(ts_list[0][j]) >= ts_threshold:
                n += 1

        Nsrc = n

        # --------------------------------- BEST FIT RA & DEC --------------------------------- !!!

        ra_list, ra_fit, dec_list, dec_fit = ([] for j in range(4))
        coord = likeObj.loadRaDec() if Ndet > 0 else None
        ra_list.append(coord[0]) if Ndet > 0 else ra_list.append([np.nan])
        dec_list.append(coord[1]) if Ndet > 0 else dec_list.append([np.nan])

        # only first elem ---!
        ra_fit.append(ra_list[0][0])
        dec_fit.append(dec_list[0][0])

        # --------------------------------- BEST FIT SPECTRAL --------------------------------- !!!

        pref_list, pref, index_list, index, pivot_list, pivot = ([] for j in range(6))
        likeObj.if_cut = if_cut
        spectral = likeObj.loadSpectral()
        index_list.append(spectral[0]) if Ndet > 0 else index_list.append([np.nan])
        pref_list.append(spectral[1]) if Ndet > 0 else pref_list.append([np.nan])
        pivot_list.append(spectral[2]) if Ndet > 0 else pivot_list.append([np.nan])

        # only first elem ---!
        index.append(index_list[0][0])
        pref.append(pref_list[0][0])
        pivot.append(pivot_list[0][0])

        # eventually cutoff ---!
        if if_cut:
            cutoff_list, cutoff = ([] for j in range(2))
            cutoff_list.append(spectral[3]) if Ndet > 0 else cutoff_list.append([np.nan])
            cutoff.append(cutoff_list[0][0])

        # --------------------------------- INTEGRATED FLUX --------------------------------- !!!

        flux_ph = []
        if Ndet > 0:
            flux_ph.append(grb.photonFluxPowerLaw(index[0], pref[0], pivot[0]))  # E (MeV)
        else:
            flux_ph.append(np.nan)

        # MISSING THE CUT-OFF OPTION ---!!!

        # --------------------------------- CLOSE LIKE XML --------------------------------- !!!

        likeObj.closeXml()

        # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

        if checks:
            print('\n!!! ---------- check trial:', count)
            print('!!! ----- check texp:', texp[i])
            print('!!! *** check Ndet:', Ndet)
            print('!!! *** check Nsrc:', Nsrc)
            print('!!! *** check ra_det:', ra_det[0])
            print('!!! *** check dec_det:', dec_det[0])
            print('!!! *** check ra_fit:', ra_fit[0])
            print('!!! *** check dec_fit:', dec_fit[0])
            print('!!! *** check flux_ph:', flux_ph[0])
            print('!!! *** check ts:', ts[0])
        
        runid = 'run0406_ID000126'
        ra = float(ra_det[0])
        dec = float(dec_det[0])
        flux = float(flux_ph[0])
        sqrt_ts = np.sqrt(float(ts[0]))
        true = SkyCoord(ra = true_coord[0]*u.deg, dec = true_coord[1]*u.deg, frame='fk5')
        dist = float(true.separation(SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='fk5')).deg)

        row = f"{runid} {count} {texp[i]} {sqrt_ts} {flux} {ra} {dec} {dist} {offset} {delay} {scaleflux} {caldb} {irf}\n"
        print(logname)
        print('!!! check row: seed %d --- texp' %i, texp[i], 's =====\n', row) if checks is True else None
        if not isfile(logname):
            hdr = 'runid seed texp sqrt_ts flux ra dec dist offset delay scaleflux caldb irf\n'
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
    if int(count) != 1:
        os.system('rm ' + p.getSimDir() + '*ebl%06d*' % count)
        os.system('rm ' + p.getSelectDir() + '*ebl%06d*' % count)
        os.system('rm ' + p.getDetDir() + '*ebl%06d*' % count)

print('\n\n!!! ================== END ================== !!!\n\n')



