# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
# *******************************************************************************

import os
import yaml

class ConfigException(Exception):
    """Alert for exception."""
    def __init__(self, message):
        super().__init__(message)

class BadConfiguration(ConfigException):
    """Alert for bad configuration."""
    def __init__(self, message):
        super().__init__(message)

class ConfigParamNotFound(ConfigException):
    """Alert when parameter is not found."""
    def __init__(self, message):
        super().__init__(message)

class Config:
    """Configure pipeline."""
    def __init__(self, cfgfile):
        configuration = open(cfgfile)
        self.cfg = yaml.load(configuration, Loader=yaml.FullLoader)
        self.cfgDesc = {
            'sections' : ['setup', 'simulation', 'analysis', 'options', 'path'],
            'setup' : ['simtype', 'runid', 'trials', 'start_count', 'scalefluxfactor'],
            'simulation' : ['caldb', 'irf', 'tobs', 'onset', 'emin', 'emax', 'roi', 'delay', 'offset', 'nruns'],
            'analysis' : ['maxsrc', 'skypix', 'skyroifrac', 'smooth', 'tool', 'type', 'blind', 'binned', 'exposure', 'usepnt', 'sgmthresh', 'cumulative', 'lightcurve', 'index'],
            'options' : ['set_ebl', 'extract_data', 'plotsky'],
            'path' : ['data', 'ebl', 'model', 'catalog']
        }
        self.validateCfg()

    def prettyd(self, d=None, indent=0):
        s = ''
        if d is None:
            d = self.cfg
        for key, value in d.items():
            s += '\t' * indent + str(key)
            if isinstance(value, dict):
                s += '\n' + self.prettyd(value, indent+1)
            else:
                s += '\t' * (indent+1) + str(value) + '\n'
        return s

    def __str__(self):
        return self.prettyd(self.cfg)
        
    def get(self, paramName):
        for sectionName in self.cfgDesc['sections']:
            if paramName in self.cfg[sectionName]:
                return self.cfg[sectionName][paramName]
        raise ConfigParamNotFound(f"Config param {paramName} not found in configuration.")

    def set(self, paramName, paramValue):
        for sectionName in self.cfgDesc['sections']:
            if paramName in self.cfg[sectionName]:
                self.cfg[sectionName][paramName] = paramValue
                return
        raise ConfigParamNotFound(f"Config param {paramName} not found in configuration.")

    def validateCfg(self):
        self.checkSections(               self.cfgDesc['sections'])
        self.checkSetupSectionParams(     self.cfgDesc['setup'])
        self.checkSimulationSectionParams(self.cfgDesc['simulation'])
        self.checkOptionsSectionParams(   self.cfgDesc['options'])
        self.checkAnalysisSectionParams(  self.cfgDesc['analysis'])
        self.checkPathSectionParams(      self.cfgDesc['path'])           
    
    def checkSections(self, sections):
        sectionMissing = set(sections) - set(self.cfg.keys())
        if len(sectionMissing) > 0:
            raise BadConfiguration(f'Configuration file sections are missing: {sectionMissing}')

    def dump(self, cfgfile_path):
        with open(cfgfile_path, 'w') as f:
            yaml.dump(self.cfg, f)

    #################
    # Setup section #
    #################

    def checkSetupSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['setup'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "setup" section are missing: {paramsMissing}')

        sectionDict = self.cfg['setup'] 
        
        # Add validations here
        simTypeValues = ['grb', 'bkg', 'skip', 'wobble', 'wilks']
        if sectionDict['simtype'] not in simTypeValues:
            raise BadConfiguration(f'simtype={sectionDict["simtype"]} is not supported. Available values: {simTypeValues}')

        # Add other validations here....

    ######################
    # Simulation section #
    ######################

    def checkSimulationSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['simulation'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "simulation" section are missing: {paramsMissing}')

        sectionDict = self.cfg['simulation'] 
        
        # Add other validations here....

    ###################
    # Options section #
    ###################

    def checkOptionsSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['options'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "options" section are missing: {paramsMissing}')

        sectionDict = self.cfg['options']

        # Add other validations here....

    ####################
    # Analysis section #
    ####################

    def checkAnalysisSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['analysis'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "analysis" section are missing: {paramsMissing}')

        sectionDict = self.cfg['analysis']

        # Add other validations here....
        toolTypeValues = ['ctools', 'gammapy', 'rtatool', 'skip']
        if sectionDict['tool'] not in toolTypeValues:
            raise BadConfiguration(f'tool={sectionDict["tool"]} is not supported. Available values: {toolTypeValues}')
        typeValues = ['1d', '3d']
        if sectionDict['type'] not in typeValues:
            raise BadConfiguration(f'type={sectionDict["type"]} is not supported. Available values: {typeValues}')

    ################
    # Path section #
    ################

    def checkPathSectionParams(self, params):
        paramsMissing = set(params) - set(self.cfg['path'])
        if len(paramsMissing) > 0:
            raise BadConfiguration(f'Configuration file params of "path" section are missing: {paramsMissing}')

        sectionDict = self.cfg['path']
        if not sectionDict['data']:
            raise BadConfiguration(f'data={sectionDict["data"]} is empty!')

        self.cfg['path']['data'] = os.path.expandvars(self.cfg['path']['data'])
        if not os.path.isdir(sectionDict['data']):
            raise BadConfiguration(f'data={sectionDict["data"]} is not a folder!')

        # Add other validations here....
        
