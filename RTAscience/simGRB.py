# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import numpy as np
import os, datetime, argparse
from pathlib import Path
from shutil import copy
from multiprocessing import Pool
from time import time

from RTAscience.lib.RTACtoolsSimulation import RTACtoolsSimulation
from RTAscience.lib.RTAManageXml import ManageXml
from RTAscience.lib.RTAUtils import get_pointing
from RTAscience.lib.RTACtoolsAnalysis import RTACtoolsAnalysis
from RTAscience.cfg.Config import Config

class GRBSimulator:

    @staticmethod
    def getConfiguration(simtype="grb", 
                               runid="run0406_ID000126", 
                               trials=1, 
                               start_count=0, 
                               irf="South_z40_average_LST_30m", 
                               tobs=18000, 
                               onset=3000, 
                               emin=0.03, 
                               emax=0.15, 
                               roi=2.5
    ):
      print("getConfiguration called!")
      configFile = "./config.tmp.yaml"
      with open(configFile, "w") as cf:
        cf.write(f"""
setup:
  simtype: {simtype}                      # grb -> src+bkg; bkg -> empty fields
  runid: {runid}           # can be all or any template (list or str) in catalog
  trials: {trials}                         # realisations per runid
  start_count: {start_count}                    # starting count for seed

simulation:
  caldb: prod3b                     # calibration database
  irf: {irf}     # istrument response function
  tobs: {tobs}                       # total obs time (s)
  onset: {onset}                        # time of bkg only a.k.a. delayed onset of burst (s)
  emin: {emin}                        # simulation minimum energy (TeV)
  emax: {emax}                        # simulation maximum energy (TeV)
  roi: {roi}                            # region of interest radius (deg)

options:
  set_ebl: True                     # uses the EBL absorbed template
  extract_data: True                # if True extracts lightcurves and spectra 

path: 
  data: $DATA                       # all data should be under this folder
  ebl: $DATA/ebl_tables/gilmore_tau_fiducial.csv
  model: $DATA/models
  catalog: $DATA/templates/grb_afterglow/GammaCatalogV1.0
"""     )
      return configFile

    def __init__(self, args):

        self.cfg = Config(args.cfgfile)

        # GRB ---!
        runid = self.cfg.get('runid')
        if type(runid) == list:
            raise ValueError('This script allows to compile only one template at a time.')

        # general ---!
        simtype = self.cfg.get('simtype')  # 'grb' -> src+bkg; 'bkg' -> empty fields

        # sim parameters ---!
        tobs = self.cfg.get('tobs')  # total obs time (s)
        onset = self.cfg.get('onset') # time of bkg only a.k.a. delayed onset of burst (s)

        datapath = self.cfg.get('data')

 
        # paths ---!
        # out_folder = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.out_folder = f"simtype_{simtype}_os_{onset}_tobs_{tobs}_irf_{self.cfg.get('irf')}_emin_{self.cfg.get('emin')}_emax_{self.cfg.get('emax')}_roi_{self.cfg.get('roi')}" 
        self.out_folder = os.path.join(datapath, 'obs', self.out_folder)

        # TODO: START_COUNT CHECK !!
        
   


        # files ---!
        ebl_table = os.path.join(datapath, 'ebl_tables/gilmore_tau_fiducial.csv')  # CSV table with EBL data
        template =  os.path.join(datapath, f'templates/{runid}.fits')  # grb FITS template data
        model_pl = os.path.join(datapath, f'models/{runid}.xml')  # grb XML template model
        tcsv = os.path.join(datapath, f'extracted_data/{runid}/time_slices.csv')  # template time bin table (to produce)
        bkg_model = os.path.join(datapath, 'models/CTAIrfBackground.xml')  # XML background model

        # Dumping the Conf object to txt file
        Path(self.out_folder).mkdir(exist_ok=True)
        dumpedConfig = os.path.join(self.out_folder, "config.yaml")


        if Path(dumpedConfig).is_file():
            Path(dumpedConfig).unlink()
        copy(args.cfgfile, str(dumpedConfig))

        # check source xml model template and create if not existing ---!
        true_coords = get_pointing(template)
        if not os.path.isfile(model_pl):
            print('Creating XML template model')
            template_pl = os.path.join(datapath, 'models/grb_file_model.xml')
            os.system('cp %s %s' % (template_pl, model_pl))
            model_xml = ManageXml(xml=model_pl)
            model_xml.setModelParameters(source=runid, parameters=('RA', 'DEC'), values=true_coords)
            del model_xml



    def start_simulation(self, cpus):
        
        trials = self.cfg.get('trials')
        sc = self.cfg.get('start_count')

        if not isinstance(sc, list):
            seeds = [i for i in range(sc, sc + trials)]
        else:
            seeds = sc
        
        print("seeds: ",seeds)

        with Pool(cpus) as p:
            times = p.map(self.simulate_trial, [ (self.cfg,seed) for seed in seeds])
        
        print(f"Times: {times}")
        if len(times) > 1:
            print(f"Mean: {np.array(times).mean()}")

    def simulate_trial(self, args_tuple):
        start_t = time()
        cfg = args_tuple[0]
        seed = args_tuple[1]
        print(f"\nSimulation for trial: {seed} started.")
        tobs = cfg.get('tobs')  # total obs time (s)
        onset = cfg.get('onset') # time of bkg only a.k.a. delayed onset of burst (s)
        tmax = tobs-onset  # total src exposure time (s)
        runid = cfg.get('runid')
        datapath = cfg.get('data')
        template =  os.path.join(datapath, f'templates/{runid}.fits')  # grb FITS template data
        pointing = get_pointing(template)
        simtype = cfg.get('simtype')  # 'grb' -> src+bkg; 'bkg' -> empty fields
        onset = cfg.get('onset') # time of bkg only a.k.a. delayed onset of burst (s)
        set_ebl = self.cfg.get('set_ebl')  # uses the EBL absorbed template

        ebl_table = os.path.join(datapath, 'ebl_tables/gilmore_tau_fiducial.csv')  # CSV table with EBL data
        template =  os.path.join(datapath, f'templates/{runid}.fits')  # grb FITS template data
        model_pl = os.path.join(datapath, f'models/{runid}.xml')  # grb XML template model
        tcsv = os.path.join(datapath, f'extracted_data/{runid}/time_slices.csv')  # template time bin table (to produce)
        bkg_model = os.path.join(datapath, 'models/CTAIrfBackground.xml')  # XML background model

        name = f'ebl{seed:06d}'
        # setup ---!
        sim = RTACtoolsSimulation()
        sim.configure(cfg)
        sim.tmax = tmax
        sim.seed = seed
        sim.nthreads = 2
        sim.pointing = pointing

        grbpath = os.path.join(datapath, 'obs', self.out_folder, runid)  # folder that will host the phlist src+bkg phlists
        bkgpath = os.path.join(datapath, 'obs', self.out_folder, 'backgrounds')  # folter that will host the bkgs only

        # check folders and create missing ones ---!
        Path(grbpath).mkdir(parents=True, exist_ok=True)
        Path(bkgpath).mkdir(parents=True, exist_ok=True)

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

            # -------------------------------------------- shift time --- !!!
            if onset != 0:
                print(f'Shift template of {onset} s after the start of the observation')
                sim.shiftTemplateTime(phlist=event_bins, time_shift=onset)

                # ------------------------------------ add background --- !!!
                print('Simulate bkg to add before the burst')
                # bkg = os.path.join(grbpath, f'{name}.fits')
                bkg = os.path.join(grbpath, f'bkg{seed:06d}.fits')            
                event_bins.insert(0, bkg)
                sim.t = [0, onset]
                sim.model = bkg_model
                sim.output = bkg
                sim.run_simulation()

            # ---------------------------- merge in single photon list ---!!!
            print('Merge in single photon list ')
            phlist = os.path.join(grbpath, f'{name}.fits')
            sim.input = event_bins
            sim.output = phlist
            sim.appendEventsSinglePhList(GTI=[0, tobs])

            sim.input = phlist
            sim.sortObsEvents()
            del sim
            print('remove template bins')
            os.system('rm ' + os.path.join(grbpath, f'{name}*tbin*'))
            if onset != 0:
                os.system('rm ' + bkg)

        # -------------------------------------------------------- BKG ---!!!
        elif simtype.lower() == 'bkg':
            print('Simulate empty fields ->',name)
            sim.seed = seed
            sim.t = [0, tobs]
            name = f'bkg{seed:06d}'
            bkg = os.path.join(bkgpath, f'{name}.fits')
            sim.model = bkg_model
            sim.output = bkg
            sim.run_simulation()

        elapsed_t = time()-start_t
        print(f"Trial {seed} took {elapsed_t} seconds")
        return (seed,elapsed_t)




if __name__=='__main__':

    parser = argparse.ArgumentParser(description='ADD SCRIPT DESCRIPTION HERE')
    parser.add_argument('--cfgfile', type=str, required=True, help="Path to the yaml configuration file")
    parser.add_argument('--ncpu', type=int, required=False, default=4, help="The required number of cpus")
    args = parser.parse_args()

    grbSim = GRBSimulator(args)
    grbSim.start_simulation(args.ncpu)
