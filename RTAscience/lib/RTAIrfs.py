# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import subprocess
import numpy as np
from os.path import join
from astropy.io import fits
from scipy.interpolate import interp1d

class RTAIrfs:
    '''
    This class allows to degrade the CTA instrument response functions.
    '''
    def __init__(self):
        # location of ctools ---!
        self.__CALDB = os.environ.get('CTOOLS')
        # files fields ---!
        self.caldb = 'prod2'  # production name in calibration database ---!
        self.irf = 'South_0.5h'  # irf ID name ---!
        # irf degradation & flux reduction ---!
        self.factor = 2


    # set CALDB var location ---!
    def setCALDB(self, path):
        '''Set path to CALDB.'''
        self.__CALDB = path
        return

    # get CALDB var value ---!
    def getCALDB(self):
        '''Get path to CALDB.'''
        return self.__CALDB

    # initialize paths for caldb degradation: directories and files ---!
    def __initCaldbIrf(self):
        '''Initialise paths and folders of nominal CALDB and a copy to degrade.'''
        nominal_irf =  f'{self.__CALDB}/share/caldb/data/cta/{self.caldb}/bcf/{self.irf}/irf_file.fits'
        degraded_irf = nominal_irf.replace('prod', 'degr')
        caldb_degr = self.caldb.replace('prod', 'degr')
        folder = f'{self.__CALDB}/share/caldb/data/cta/'
        nominal_cal = join(folder, self.caldb)
        degraded_cal = join(folder, caldb_degr)
        return folder, nominal_cal, nominal_irf, degraded_cal, degraded_irf

    # updates the degraded caldb index by replacing all "prod" references with "degr" ---!
    def __updateCaldbIndex(self, index):
        '''Updates the CALDB index.'''
        # read content ---!
        with open(index, 'r', encoding="ISO-8859-1") as f:
            filedata = f.read()
        # Replace the target keyword ---!
        filedata = filedata.replace('prod', 'degr').replace('PROD', 'DEGR')
        # Write the file out again ---!
        with open(index, 'w', encoding="ISO-8859-1") as f:
            f.write(filedata)
        return

    # create copy of caldb and corresponding caldb.inx file ---!
    def __mockNominalCaldb(self, nominal_cal, nominal_irf, degraded_cal, degraded_irf):
        '''Generates a copy of the nominal CALDB.'''
        if not os.path.isdir(degraded_cal):
            os.mkdir(degraded_cal)
        if not os.path.isfile(join(degraded_cal,'caldb.indx')):
            os.system(f"cp {join(nominal_cal, 'caldb.indx')} {join(degraded_cal, 'caldb.indx')}")
            # update caldb.indx file ---!
            self.__updateCaldbIndex(join(degraded_cal, 'caldb.indx'))
        if not os.path.isdir(join(degraded_cal, 'bcf')):
            os.mkdir(join(degraded_cal, 'bcf'))
        if not os.path.isdir(join(degraded_cal, 'bcf', self.irf)):
            os.mkdir(join(degraded_cal, 'bcf', self.irf))
        if os.path.isfile(degraded_irf):
            os.system(f'rm {degraded_irf}')
        if not os.path.isfile(degraded_irf):
            os.system(f'cp {nominal_irf} {degraded_irf}')
        return

    # change permission to 777 and ask for password if user id not in idlist param ---!
    def __openPermission(self, path, idlist=(0,1126,1001)):
        '''Grants writing permission to the CALDB folder.'''
        if os.geteuid() in idlist:
            subprocess.run(['chmod', '-R', '777', path], check=True)
        else:
            subprocess.run(['sudo', 'chmod', '-R', '777', path], check=True)
        return

    # change permission to 755 and ask for password if user id not in idlist param ---!
    def __closePermission(self, path, idlist=(0,1126)):
        '''Removes writing permission to the CALDB folder'''
        if os.geteuid() in idlist:
            subprocess.run(['chmod', '-R', '755', path], check=True)
        else:
            subprocess.run(['sudo', 'chmod', '-R', '755', path], check=True)
        return

    # degrade Aff by self.factor (for now only scalar is implemented) ---!
    def __degradeAeff(self, nominal_irf, degraded_irf, r=False):
        '''Modifies the AEFF matrix by a factor (scalar).'''
        # initialise ---!
        inv = 1 / self.factor
        extension = 'EFFECTIVE AREA'
        field = 4
        with fits.open(nominal_irf) as hdul:
            elo = np.array(hdul[extension].data.field(0)[:].astype(float)[0])
            ehi = np.array(hdul[extension].data.field(1)[:].astype(float)[0])
            e = elo + 0.5*(ehi - elo)
            tlo = np.array(hdul[extension].data.field(2)[:].astype(float)[0])
            thi = np.array(hdul[extension].data.field(3)[:].astype(float)[0])
            theta = tlo + 0.5*(thi - tlo)
            aeff = np.array(hdul[extension].data.field(field)[:].astype(float)[0])
        # effective area multiplied by inv of factor ---!
        a = np.where(np.array([i * inv for i in aeff]) is np.nan, 0., np.array([i * inv for i in aeff]))
        # degrade and save new ---!
        with fits.open(degraded_irf, mode='update') as hdul:
            hdul[extension].data.field(field)[:] = a
            # save changes ---!
            hdul.flush()
        # return only if bkg counts must be degraded ---!
        if not r:
            return
        else:
            return aeff, a, theta, e

    # degrade bkg counts by normalise for aeff nominal and multiply times aeff degraded ---!
    def __degradeBkg(self, nominal_irf, degraded_irf, aeff=True):
        '''Modifies the BKG matrix by a factor (scalar).'''
        # degrade Aeff (only if True) and get its returns ---!
        if not aeff:
            tmp = self.factor
            self.factor = 1
        aeff_nom, aeff_deg, theta, e_aeff = self.__degradeAeff(nominal_irf=nominal_irf, degraded_irf=degraded_irf, r=True)
        if not aeff:
            self.factor = tmp
        # initialise ---!
        extension = 'BACKGROUND'
        field = 6
        with fits.open(nominal_irf) as hdul:
            xlo = np.array(hdul[extension].data.field(0)[:].astype(float)[0])
            xhi = np.array(hdul[extension].data.field(1)[:].astype(float)[0])
            x = xlo + 0.5*(xhi - xlo)
            ylo = np.array(hdul[extension].data.field(2)[:].astype(float)[0])
            yhi = np.array(hdul[extension].data.field(3)[:].astype(float)[0])
            y = ylo + 0.5*(yhi - ylo)
            elo = np.array(hdul[extension].data.field(4)[:].astype(float)[0])
            ehi = np.array(hdul[extension].data.field(5)[:].astype(float)[0])
            e_bkg = elo + 0.5*(ehi - elo)
            bkg = np.array(hdul[extension].data.field(field)[:].astype(float)[0])
        # spatial pixel/deg conversion factor ---!
        conv_factor = (xhi.max() - xlo.min()) / theta.max()
        # interpolated Aeff via energy grid ---!
        nominal_interp, degraded_interp = ([[]*i for i in range(len(theta))] for i in range(2))
        for i in range(len(theta)):
            fnom = interp1d(e_aeff[:], aeff_nom[i,:])
            nominal_interp[i].append(fnom(e_bkg[:]))
            fdeg = interp1d(e_aeff[:], aeff_deg[i,:])
            degraded_interp[i].append(fdeg(e_bkg[:]))
        # flatten list of theta interpolations (theta array of energy frames) ---!
        nominal_interp = np.array([item for sublist in nominal_interp for item in sublist])
        degraded_interp = np.array([item for sublist in degraded_interp for item in sublist])
        # empty copy of bkg tensor ---!
        b = np.empty_like(bkg)
        for idf, frame in enumerate(bkg[:,0,0]):
            for idx, xpix in enumerate(bkg[idf,:,0]):
                for idy, ypix in enumerate(bkg[idf,idx,:]):
                    # find radius in degrees ---!
                    r = np.sqrt((0 - xpix)**2 + (0 - ypix)**2)
                    rdegree = r * conv_factor
                    # find corresponding theta index ---!
                    angle = min(theta, key=lambda x:abs(x-rdegree))
                    idtheta = np.where(np.isin(theta[:], angle))
                    # degrade the background count for frame/x/y point ---!
                    if nominal_interp[idtheta,idf] == 0.:
                        b[idf, idx, idy] = 0.
                    else:
                        b[idf,idx,idy] = bkg[idf,idx,idy] / nominal_interp[idtheta,idf] * degraded_interp[idtheta,idf]
        # save to new ---!
        with fits.open(degraded_irf, mode='update') as hdul:
            hdul[extension].data.field(field)[:] = b
            # save changes ---!
            hdul.flush()
        return

    # degrade IRFs via Effective Area and/or Background ---!
    def degradeIrf(self, bkg=True, aeff=True, mod_permission=False):
        '''From a nominal CALDB generates a degraded copy.'''
        # initialize ---!
        folder, nominal_cal, nominal_irf, degraded_cal, degraded_irf = self.__initCaldbIrf()
        # open all folder permission ---!
        if mod_permission:
            self.__openPermission(path=folder)
        # create degr caldb path if not existing ---!
        self.__mockNominalCaldb(nominal_cal=nominal_cal, nominal_irf=nominal_irf, degraded_cal=degraded_cal, degraded_irf=degraded_irf)
        # close all folder permission and open only degraded caldb permission ---!
        if mod_permission:
            self.__closePermission(path=folder)
            self.__openPermission(path=degraded_cal)
        # degradation aeff ---!
        if not bkg:
            self.__degradeAeff(nominal_irf=nominal_irf, degraded_irf=degraded_irf)
        # degradation bkg counts ---!
        else:
            self.__degradeBkg(nominal_irf=nominal_irf, degraded_irf=degraded_irf, aeff=aeff)
        # close degraded caldb permission ---!
        if mod_permission:
            self.__closePermission(degraded_cal)
        # update caldb ---!
        self.caldb = self.caldb.replace('prod', 'degr')
        return