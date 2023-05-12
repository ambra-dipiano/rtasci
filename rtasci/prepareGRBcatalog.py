# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import argparse
from os.path import isdir, join
from rtasci.lib.RTACtoolsSimulation import RTACtoolsSimulation
from rtasci.lib.RTAUtils import get_pointing
from rtasci.lib.RTAManageXml import ManageXml
from rtasci.cfg.Config import Config

parser = argparse.ArgumentParser(description='This script extracts spectra and lightcurves from the GRB templates, in order to prepare all required files for the simulation.')
parser.add_argument('-f', '--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
args = parser.parse_args()

cfg = Config(args.cfgfile)

# GRB ---!
if cfg.get('runid') == 'all':
    runids = [f.replace('.fits', '') for f in os.listdir(cfg.get('catalog')) if '_ebl' not in f and os.path.isfile(join(cfg.get('catalog'), f))]
elif type(cfg.get('runid')) == str:
    runids = [cfg.get('runid')]
else:
    runids = cfg.get('runid')
runids = sorted(runids)

# check scalefluxfactor ---!
if cfg.get('scalefluxfactor') == None:
    raise ValueError('The parameter "scalefluxfactor" must be int or float. If you want to use the nominal template, please set scalefluxfactor=1.')

# conditions control ---!
set_ebl = cfg.get('set_ebl')  # uses the EBL absorbed template
# paths ---!
if '$' in cfg.get('data'):
    datapath = os.path.expandvars(cfg.get('data'))
else:
    datapath = cfg.get('data')  
if not isdir(datapath):
    raise ValueError('Please specify a valid path')
if not isdir(join(datapath, 'obs')):
    os.mkdir(join(datapath, 'obs'))

# global files ---!
ebl_table = os.path.expandvars(cfg.get('ebl')) 
pl_template = join(os.path.expandvars(cfg.get('model')), 'grb_file_model.xml' )
if not os.path.isfile(pl_template):
    raise ValueError(f'PL template {pl_template} not found')
if not os.path.isfile(ebl_table):
    raise ValueError(f'EBL table {ebl_table} not found')


# ------------------------------------------------ loop over runids

for runid in runids:
    print(f'Processing runid: {runid}')
    # grb path ---!
    grbpath = join(datapath, 'obs', runid)  # folder that will host the phlist 
    if not isdir(grbpath):
        os.mkdir(grbpath)
    
    # grb files ---!
    template =  join(os.path.expandvars(cfg.get('catalog')), f'{runid}.fits')  # grb FITS template data
    model_pl = pl_template.replace('grb_file_model.xml', f'{runid}.xml')  # grb XML template model
    if not isdir(join(datapath, f'extracted_data')):
        os.mkdir(join(datapath, f'extracted_data'))
    tcsv = join(datapath, f'extracted_data/{runid}/time_slices.csv')  # grb template time grid
    if not os.path.isfile(template):
        raise ValueError(f'Template {runid} FITS not found')
    if not isdir(join(datapath, f'extracted_data/{runid}')):
        os.mkdir(join(datapath, f'extracted_data/{runid}'))
    if not os.path.isfile(model_pl):
        print(f'Creating {runid} XML model')
        os.system(f'cp {pl_template} {model_pl}')
        os.system(f'chmod +wr {model_pl}')
        model_xml = ManageXml(xml=model_pl)
        model_xml.setModelParameters(source=runid, parameters=('RA', 'DEC'), values=get_pointing(template))
        del model_xml

    sim = RTACtoolsSimulation()
    sim.template = template
    sim.model = model_pl
    # add EBL to template ---!
    if set_ebl:
        sim.table = ebl_table  
        sim.zfetch = True
        if not sim.checkEBLinFITS():
            print('Computing EBL absorption')
            sim.set_ebl = False
            sim.addEBLtoFITS(template, ext_name='EBL-ABS. SPECTRA')
        sim.template = template
    sim.set_ebl = set_ebl
    # load template ---!
    if cfg.get('extract_data'):
        sim.extract_spectrum = True
        print('Creating lightcurves and spectra')
    sim.loadTemplate(source_name=runid, return_bin=False, data_path=join(datapath, f'extracted_data/{runid}'), scalefluxfactor=cfg.get('scalefluxfactor'))

