# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import sys
import os
import yaml
import numpy as np
from lib.RTACtoolsSimulation import RTACtoolsSimulation
from lib.RTAUtils import get_pointing
from lib.RTAManageXml import ManageXml

cfgfile = sys.argv[1]
pypath = str(os.path.dirname(os.path.abspath(__file__)))  
configuration = open(os.path.join(pypath, cfgfile) )
cfg = yaml.load(configuration, Loader=yaml.FullLoader)

# GRB ---!
if cfg['setup']['runid'] == 'all':
    runids = [f.replace('.fits', '') for f in os.listdir(cfg['path']['catalog']) if '_ebl' not in f and os.path.isfile(os.path.join(cfg['path']['catalog'], f))]
elif type(cfg['setup']['runid']) == str:
    runids = [cfg['setup']['runid']]
else:
    runids = cfg['setup']['runid']

# conditions control ---!
set_ebl = True  # uses the EBL absorbed template
# paths ---!
if '$' in cfg['path']['data']:
    datapath = os.path.expandvars(cfg['path']['data'])
else:
    datapath = cfg['path']['data']  
# global files ---!
ebl_table = os.path.join(datapath, 'ebl_tables/gilmore_tau_fiducial.csv') 
pl_template = os.path.join(datapath, 'models/grb_file_model.xml')
if not os.path.isfile(pl_template):
    raise ValueError(f'PL template {pl_template} not found')
if not os.path.isfile(ebl_table):
    raise ValueError(f'EBL table {ebl_table} not found')


# ------------------------------------------------ loop over runids

for runid in runids:
    print(f'Processing runid: {runid}')
    # grb path ---!
    grbpath = os.path.join(datapath, 'obs', runid)  # folder that will host the phlist 
    if not os.path.isdir(datapath):
        raise ValueError('Please specify a valid path')
    if not os.path.isdir(grbpath):
        if not os.path.isdir(os.path.join(datapath, 'obs')):
            os.mkdir(os.path.join(datapath, 'obs'))
        os.mkdir(grbpath)
    # grb files ---!
    template =  os.path.join(datapath, f'templates/grb_afterglow/GammaCatalogV1.0/{runid}.fits')  # grb FITS template data
    model_pl = os.path.join(datapath, f'models/{runid}.xml')  # grb XML template model
    tcsv = os.path.join(datapath, f'extracted_data/{runid}/time_slices.csv')  # grb template time grid
    if not os.path.isfile(template):
        raise ValueError(f'Template {runid} FITS not found')
    if not os.path.isdir(os.path.join(datapath, f'extracted_data/{runid}')):
        os.mkdir(os.path.join(datapath, f'extracted_data/{runid}'))
    true_coords = get_pointing(template)
    if not os.path.isfile(model_pl):
        print(f'Creating {runid} XML model')
        os.system(f'cp {pl_template} {model_pl}')
        os.system(f'chmod +wr {model_pl}')
        model_xml = ManageXml(xml=model_pl)
        model_xml.setModelParameters(source=runid, parameters=('RA', 'DEC'), values=true_coords)
        del model_xml

    sim = RTACtoolsSimulation()
    sim.template = template
    sim.model = model_pl
    # add EBL to template ---!
    if set_ebl:
        sim.table = ebl_table  
        sim.zfetch = True
        sim.set_ebl = False
        if not os.path.isfile(template.replace('.fits', '_ebl.fits')):
            print('Computing EBL absorption')
            sim.addEBLtoFITS(template.replace('.fits', '_ebl.fits'), ext_name='EBL-ABS. SPECTRA')
        sim.set_ebl = set_ebl
        sim.template = template.replace('.fits', '_ebl.fits')
    # load template ---!
    if not os.path.isfile(tcsv):
        sim.extract_spectrum = True
        print('Creating lightcurves and spectra')
    tbin_stop = sim.loadTemplate(source_name=runid, return_bin=True, data_path=os.path.join(datapath, f'extracted_data/{runid}'))