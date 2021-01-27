# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import yaml
import numpy as np
from lib.RTACtoolsSimulation import RTACtoolsSimulation

cfgfile = sys.argv[1]
pypath = str(os.path.dirname(os.path.abspath(__file__)))  
configuration = open(join(pypath, cfgfile) )
cfg = yaml.load(configuration, Loader=yaml.FullLoader)

# files ---!
ebl_table = os.path.join(cfg['path']['data'], 'ebl_tables/gilmore_tau_fiducial.csv')  # CSV table with EBL data
template =  os.path.join(cfg['path']['data'], f'templates/{runid}.fits')  # grb FITS template data
model_pl = os.path.join(cfg['path']['data'], f'models/{runid}.xml')  # grb XML template model
tcsv = os.path.join(cfg['path']['data'], f'extracted_data/{runid}/time_slices.csv')  # template time bin table (to produce)
bkg_model = os.path.join(cfg['path']['data'], 'models/CTAIrfBackground.xml')  # XML background model


# check destination folder
if not os.path.isdir(os.path.join(cfg['path']['data'], f'extracted_data/{runid}')):
    os.mkdir(os.path.join(cfg['path']['data'], f'extracted_data/{runid}'))



# check source xml model template and create if not existing ---!
if not os.path.isfile(model_pl):
    print('Creating XML template model')
    template_pl = os.path.join(cfg['path']['data'], 'models/grb_file_model.xml')
    os.system('cp %s %s' % (template_pl, model_pl))
    model_xml = ManageXml(xml=model_pl)
    model_xml.setModelParameters(source=runid, parameters=('RA', 'DEC'), values=true_coord)
    del model_xml



sim = RTACtoolsSimulation()
sim.seed = count
sim.nthreads = nthreads
sim.pointing = pointing
sim.roi = roi
sim.e = [emin, emax]
sim.tmax = tmax
sim.caldb = caldb
sim.irf = irf

# -------------------------------------------------------- GRB ---!!!
if simtype.lower() == 'grb':
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
    sim.extract_spectrum = extract_spectrum
    tbin_stop = sim.loadTemplate(source_name=runid, return_bin=True, data_path=os.path.join(cfg['path']['data'], f'extracted_data/{runid}'))