# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
# *******************************************************************************

from rtasci.cfg.Config import Config

class RTACtoolsBase:

    def __init__(self):
        self.caldb = 'prod2'  # production name in calibration database ---!
        self.irf = 'South_0.5h'  # irf ID name ---!
        # data ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.roi = 5  # region of indeterest (deg) ---!
    
    def configure(self, cfg: Config):
        '''Sets common parameters for analysis and simulation: caldb, irf, energy range, field of view.'''
        self.caldb = cfg.get('caldb')
        self.irf = cfg.get('irf')
        self.e = [cfg.get('emin'), cfg.get('emax')] 
        self.roi = cfg.get('roi')