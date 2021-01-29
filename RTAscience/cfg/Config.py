# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
# *******************************************************************************

import yaml

class ConfigException(Exception):
    def __init__(self, message):
        super().__init__(message)

class BadConfiguration(ConfigException):
    def __init__(self, message):
        super().__init__(message)

class ConfigParamNotFound(ConfigException):
    def __init__(self, message):
        super().__init__(message)

class Config:
    def __init__(self, cfgfile):
        configuration = open(cfgfile)
        self.cfg = yaml.load(configuration, Loader=yaml.FullLoader)
        self.cfgDesc = {
            'sections' : ['setup', 'simulation', 'options', 'path'],
            'setup' : ['simtype', 'runid', 'trials', 'start_count'],
            'simulation' : ['caldb', 'irf', 'tobs', 'emin', 'emax', 'roi'],
            'options' : ['set_ebl', 'extract_data'],
            'path' : ['data', 'ebl', 'model', 'catalog']
        }
        self.validateCfg()

    def get(self, paramName):
        for sectionName in self.cfgDesc['sections']:
            if paramName in self.cfg[sectionName]:
                return self.cfg[sectionName][paramName]
        raise ConfigParamNotFound(f"Config param {paramName} not found in configuration.")

    def validateCfg(self):
        self.checkSections(               self.cfgDesc['sections'])
        self.checkSetupSectionParams(     self.cfgDesc['setup'])
        self.checkSimulationSectionParams(self.cfgDesc['simulation'])
        self.checkOptionsSectionParams(   self.cfgDesc['options'])
        self.checkPathSectionParams(      self.cfgDesc['path'])           
    
    def checkSections(self, sections):
        sectionMissing = set(sections) - set(self.cfg.keys())
        if len(sectionMissing) > 0:
            raise BadConfiguration(f'Configuration file sections are missing: {sectionMissing}')

    #################
    # Setup section #
    #################

    def checkSetupSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['setup'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "setup" section are missing: {paramsMissing}')

        sectionDict = self.cfg['setup'] 
        
        # Add validations here
        
        simTypeValues = ['grb', 'bkg']
        if sectionDict['simtype'] not in simTypeValues:
            raise BadConfiguration(f'simtype={sectionDict["simtype"]} is not supported. Available values: {simTypeValues}')

        # Add other validations here....
        # ....

    ######################
    # Simulation section #
    ######################

    def checkSimulationSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['simulation'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "simulation" section are missing: {paramsMissing}')

        sectionDict = self.cfg['simulation'] 
        

        # Add other validations here....
        # ....        


    ###################
    # Options section #
    ###################

    def checkOptionsSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['options'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "simulation" section are missing: {paramsMissing}')

        sectionDict = self.cfg['options']

        # Add other validations here....
        # ....        



    ################
    # Path section #
    ################

    def checkPathSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['path'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "path" section are missing: {paramsMissing}')

        sectionDict = self.cfg['path']

        # Add other validations here....
        
