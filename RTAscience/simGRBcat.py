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
from lib.RTAManageXml import ManageXml
from lib.RTAUtils import get_pointing

# GRB ---!
runid = ['run0406_ID000126', 'run0852_ID000211']
# general ---!
trials = 2  # trials
count = 0  # starting count
nthreads = 2
# sim parameters ---!
caldb = 'prod3b'  # calibration database
irf = 'South_z40_average_LST_30m'  # istrument response function
tmax = 12e2  # total obs time (s)
emin = 3e-2  # simulation minimum energy (TeV)
emax = 1.  # simulation maximum energy (TeV)
roi = 2.5  # region of interest radius (deg)
# conditions control ---!
set_ebl = True  # uses the EBL absorbed template
# paths ---!
pypath = str(os.path.dirname(os.path.abspath(__file__)))  
datapath = pypath.replace('cta-sag-sci', 'DATA')  # all data should be under this folder
grbpath = os.path.join(datapath, 'obs', runid)  # folder that will host the phlist 
# files ---!
ebl_table = os.path.join(datapath, 'ebl_tables/gilmore_tau_fiducial.csv')  # CSV table with EBL data
template =  os.path.join(datapath, f'templates/{runid}.fits')  # grb FITS template data
model_pl = os.path.join(datapath, f'models/{runid}.xml')  # grb XML template model
tcsv = os.path.join(datapath, f'extracted_data/{runid}/time_slices.csv')  # times table 
bkg_model = os.path.join(datapath, 'models/CTAIrfBackground.xml')  # XML background model
# check folders and create missing ones ---!
if not os.path.isdir(grbpath):
    if not os.path.isdir(grbpath.replace(runid, '')):
        os.mkdir(grbpath.replace(runid, ''))
    os.mkdir(grbpath)
if not os.path.isdir(os.path.join(datapath, f'extracted_data/{runid}')):
    os.mkdir(os.path.join(datapath, f'extracted_data/{runid}'))

# check source xml model template and create if not existing ---!
if not os.path.isfile(model_pl):
    print('Creating XML template model')
    template_pl = os.path.join(datapath, 'models/grb_file_model.xml')
    os.system('cp %s %s' % (template_pl, model_pl))
    model_xml = ManageXml(xml=model_pl)
    model_xml.setModelParameters(source=runid, parameters=('RA', 'DEC'), values=true_coord)
    del model_xml

# ---------------------------------------------------- trials --- !!!
pointing = get_pointing(template)
while count < trials:
    count += 1
    name = f'ebl{count:06d}'
    # setup ---!
    sim = RTACtoolsSimulation()
    sim.seed = count
    sim.nthreads = nthreads
    sim.pointing = pointing
    sim.roi = roi
    sim.e = [emin, emax]
    sim.tmax = tmax
    sim.caldb = caldb
    sim.irf = irf

    print(f'Simulate GRB + BKG with onset = {onset} s')
    sim.template = template
    sim.model = model_pl
    # add EBL to template ---!
    if set_ebl:
        print('Computing EBL absorption')
        sim.table = ebl_table  
        sim.zfetch = True
        sim.set_ebl = False
        if not os.path.isfile(template.replace('.fits', '_ebl.fits')):
            sim.addEBLtoFITS(template.replace('.fits', '_ebl.fits'), ext_name='EBL-ABS. SPECTRA')
        sim.set_ebl = set_ebl
        sim.template = template.replace('.fits', '_ebl.fits')
    # load template ---!
    if not os.path.isfile(tcsv):
        sim.extract_spectrum = True
    tbin_stop = sim.loadTemplate(source_name=runid, return_bin=True, data_path=os.path.join(datapath, f'extracted_data/{runid}'))

    event_bins = []
    # get time grid ---!
    sim.table = tcsv
    tgrid = sim.getTimeSlices(GTI=(0, tmax)) 
    # ----------------------------------------------- simulate ---!!!
    for i in range(tbin_stop):
        sim.t = [tgrid[i], tgrid[i + 1]]
        sim.model = os.path.join(datapath, f'extracted_data/{runid}/{runid}_tbin{i:02d}.xml')
        event = os.path.join(grbpath, f'{name}_tbin{i:02d}.fits')
        event_bins.append(event)
        sim.output = event
        sim.run_simulation()

    # ---------------------------- merge in single photon list ---!!!
    phlist = os.path.join(grbpath, f'{name}.fits')
    sim.input = event_bins
    sim.output = phlist
    sim.appendEventsSinglePhList(GTI=[0, tmax])

    del sim
    print('os.remove template bins')
    os.system('rm ' + os.path.join(grbpath, f'{name}*tbin*'))
