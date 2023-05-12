# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import gammalib
import ctools
import os.path
import csv
import re
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.interpolate import interp1d

# create observation list with gammalib ---!
def make_obslist(obslist, items, names, instruments='CTA'):
    '''Generates an observation list XML file using gammalib.'''
    if type(items) != type(list()):
        items = [items]
    if type(names) != type(list()):
        names = [names for i in range(len(items))]
    if type(instruments) != type(list()):
        instruments = [instruments for i in range(len(items))]
    xml = gammalib.GXml()
    obslib = xml.append('observation_list title="observation library"')
    for i, item in enumerate(items):
        obs = obslib.append(f'observation name="{names[i]}" id="{i+1:02d}" instrument="{instruments[i]}"')
        obs.append(f'parameter name="EventList" file="{item}"')
    xml.save(obslist)
    del xml
    return 

class RTACtoolsSimulation():
    '''
    This class allows to: 1) compute the EBL absorption from a csv data table and add it to the template; 2) extract spectra, lightcuves and time slices from the template (the flux values can also be normalised by a factor); 3) merge bins of the template simulation in a single photon list; 4) perform simulations using ctoobssim from ctools software package.
    '''
    def __init__(self):
        # files fields ---!
        self.model, self.template, self.table = (str() for i in range(3))
        self.output, self.input = (str() for i in range(2))
        self.caldb = 'prod2'  # caldb (str) ---!
        self.irf = 'South_0.5h'  # irf (Str) ---!
        # condition control ---!
        self.set_ebl = True  # set/unset EBL absorption feature ---!
        self.extract_spectrum = False  # set/unset spectra extraction feature ---!
        self.plot = False  # option for retrieving plotting values ---!
        self.zfetch = False  # set/unset automatic fetching of redshift ---!
        self.set_debug = False  # set/unset debug mode for ctools ---!
        self.set_log = True  # set/unset logfiles for ctools ---!
        # data ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.fov = 5  # region of interest (deg) ---!
        self.tmax = 1800  # maximum exposure time needed (s) ---!
        self.t = [0, 1800]  # time range (s/MJD) ---!
        self.pointing = [83.63, 22.01]  # RA/DEC or GLON/GLAT (deg) ---!
        # ctools miscellaneous ---!
        self.edisp = False  # set/unset edisp
        self.seed = 1  # MC seed ---!
        self.nthreads = 1  # run in parallel 
        # ebl specifics ---!
        self.z = 0.1  # redshift value ---!
        self.z_ind = 1  # redshift value index ---!
        # fits extension array ---!
        self.__time, self.__energy, self.__spectra, self.__ebl = (float() for i in range(4))

    # open and close the FITS files ---!
    def __openFITS(self):
        '''Opens FITS file.'''
        hdul = fits.open(self.template)
        return hdul
    def __closeFITS(self, hdul):
        '''Closes FITS file.'''
        hdul.close()
        return

    # retrive FITS data ---!
    def __getFitsData(self):
        '''Loads time, energy, spectra and (if present) absorbed spectra from a FITS file.'''
        hdul = self.__openFITS()
        self.__energy = np.array(hdul[1].data)
        self.__time = np.array(hdul[2].data)
        self.__Nt = len(self.__time)
        self.__Ne = len(self.__energy)
        self.__spectra = np.array(hdul[3].data)
        if self.set_ebl:
            try:
                self.__ebl = np.array(hdul[4].data)
            except:
                raise IndexError('Template extensions out of range. Unable to load EBL absorbed spectra.')
        self.__closeFITS(hdul)
        return

    # check if EBL extention already in template ---!
    def checkEBLinFITS(self, ext_name='EBL-ABS. SPECTRA'):
        '''Checks if specified extension is present in FITS file.'''
        hdul = self.__openFITS()
        try:
            ext = hdul[ext_name]
            return True
        except KeyError:
            return False

    # load csv table in pandas DataFrame and drop NaN values---!
    def __openCSV(self):
        '''Opens a CSV data file.'''
        df = pd.read_csv(self.table)
        df.dropna()
        return df

    # retrive csv data ---!
    def __getEBLfromCSV(self):
        '''Gets optical depth values from a CSV table.'''
        df = self.__openCSV()
        cols = list(df.columns)
        tau_table = np.array(df[cols[self.z_ind]])
        E = np.array(df[cols[0]]) / 1e3  # MeV --> GeV ---!
        return tau_table, E

    # retrive csv temporal bin grid of the template in use and return the necessary slice ---!
    def getTimeSlices(self, GTI, return_bins=False):
        '''Gets the time slices from a GRB afterglow template, within a given interval.'''
        self.__getFitsData()
        df = self.__openCSV()
        cols = list(df.columns)
        self.__time = np.append(0, np.array(df[cols[1]]))
        #self.__time = np.array(df[cols[1]])
        bin_start = 0
        bin_stop = 1
        for i in range(len(self.__time)):
            if self.__time[i] < GTI[0]:
                bin_start += 1
                continue
            elif self.__time[i] >= GTI[1]:
                self.__time[i] = GTI[1]
                bin_stop += i
                break
        if bin_stop <= self.__Nt:
            time_slice = slice(bin_start, bin_stop + 1)
        else:
            time_slice = slice(bin_start, bin_stop)
        if not time_slice:
            raise ValueError('Invalid GTI: cannot extract time slices')
        tgrid = self.__time[time_slice]
        tgrid[0] = GTI[0]
        if not return_bins:
            return tgrid
        else:
            return tgrid, bin_start, bin_stop

    # compute the EBL absorption ---!
    def __addEBL(self, unit='MeV'):
        '''Computes the EBL absorption.'''
        self.__getFitsData()
        tau_table, E = self.__getEBLfromCSV()
        if unit == 'GeV':
            E *= 1e3
        elif unit == 'TeV':
            E *= 1e6
        # interpolate linearly handling NaNs/inf/zeroes ---!
        with np.errstate(invalid='raise'):
            interp = interp1d(E, tau_table, bounds_error=False)
        tau = np.array(interp(self.__energy))
        self.__ebl = np.empty_like(self.__spectra)
        # compute absorption ---!
        for i in range(self.__Nt):
            for j in range(self.__Ne):
                self.__ebl[i][j] = self.__spectra[i][j] * np.exp(-tau[j])
        # if required return values to plot ---!
        if self.plot:
            return E, tau_table, self.__energy, tau
        else:
            return

    # retrive redshift, find nearest column then access its index ---!
    def __zfetch(self):
        '''Retrives the optical depth values from a table, according to redshift.'''
        hdul = self.__openFITS()
        # fetch z from the template and chose the table column with min distance from it 
        z = hdul[0].header['REDSHIFT']
        with open(self.table, 'r') as f:
            reader = csv.reader(f)
            hdr = next(reader)
        zlist = []
        # load only the redshift columns 
        for el in hdr:
            zlist.append(re.sub('[^0-9,.]', '', el))
        zlist.remove('')
        zlist = [float(i) for i in zlist]
        # find nearest ---!
        self.z = min(zlist, key=lambda x:abs(x-z))
        self.z_ind = zlist.index(self.z) +1
        return

    # add EBL extension to a FITS template ---!
    def addEBLtoFITS(self, template_ebl, ext_name='EBL ABS. SPECTRA', unit='MeV'):
        '''Adds the EBL absorbed spectra to the tempalte.'''
        hdul = self.__openFITS()
        if self.zfetch:
            self.__zfetch()
        # if required retrive values to plot ---!
        if self.plot:
            x, y, x2, y2 = self.__addEBL(unit=unit)
        else:
            self.__addEBL(unit=unit)
        # update fits ---!
        hdu = fits.BinTableHDU(name=ext_name, data=self.__ebl)
        header = hdu.header
        header.set('UNITS', 'ph/cm2/s/GeV', ' ')
        hdu = fits.BinTableHDU(name=ext_name, data=self.__ebl, header=header)
        hdul.append(hdu)
        # save to new ---!
        if os.path.isfile(template_ebl):
            os.remove(template_ebl)
        hdul.writeto(template_ebl, overwrite=True)
        self.__closeFITS(hdul)
        # if required return values to plot ---!
        if self.plot:
            return x, y, x2, y2
        else:
            return

    # extract template spectra, create xml model files and time slices csv file ---!
    def __extractSpectrumAndModelXML(self, source_name, time_slice_name='time_slices.csv', data_path=None, scalefluxfactor=1):
        '''Generates spectra, lightcurves and time slices of a template.'''
        # time slices table ---!
        if data_path is None:
            raise ValueError('please specify a valid path')
        table = os.path.join(data_path, time_slice_name)
        if os.path.isfile(table):
            os.remove(table)
        with open(table, 'w+') as tab:
            tab.write('#bin,tmax_bin')
        # spectra and models ---!
        for i in range(self.__Nt):
            filename = os.path.join(data_path, f'spec_tbin{i:02d}.out')
            if os.path.isfile(filename):
                os.remove(filename)
            # time slices table ---!
            with open(table, 'a') as tab:
                tab.write('\n' + str(i) + ', ' + str(self.__time[i][0]))
            # spectra ---!
            with open(filename, 'a+') as f:
                for j in range(self.__Ne):
                    # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
                    if self.set_ebl:
                        f.write(str(self.__energy[j][0] * 1000.0) + ' ' + str(self.__ebl[i][j] / 1000.0 / scalefluxfactor) + "\n")
                    else:
                        f.write(str(self.__energy[j][0] * 1000.0) + ' ' + str(self.__spectra[i][j] / 1000.0 / scalefluxfactor) + "\n")
            # xml models ---!
            os.system('cp ' + str(self.model) + ' ' + str(os.path.join(data_path, f'{source_name}_tbin{i:02d}.xml')))
            s = open(os.path.join(data_path, f'{source_name}_tbin{i:02d}.xml')).read()
            s = s.replace('data/spec', f'spec_tbin{i:02d}')
            with open(os.path.join(data_path, f'{source_name}_tbin{i:02d}.xml'), 'w') as f:
                f.write(s)
        return

    # read template and return tbin_stop containing necessary exposure time coverage ---!
    def loadTemplate(self, source_name, return_bin=False, data_path=None, scalefluxfactor=1):
        '''Loads template data (spectra, lightcurves and time slices).'''
        self.__getFitsData()
        # time grid ---!
        t = [0.0 for x in range(self.__Nt + 1)]
        for i in range(self.__Nt - 1):
            t[i + 1] = self.__time[i][0] + (self.__time[i + 1][0] - self.__time[i][0]) / 2
        # tmax in last bin ---!
        t[self.__Nt] = self.__time[self.__Nt - 1][0] + (self.__time[self.__Nt - 1][0] - t[self.__Nt - 1])
        # stop the second after higher tmax ---!
        if self.tmax != None:
            tbin_stop = 1
            for bin in range(len(t)):
                if t[bin] <= self.tmax:
                    tbin_stop += 1
                else:
                    continue
        else:
            raise ValueError('Total exposure time longer than template temporal evolution.')
        # energy grid ---!
        en = [1.0 for x in range(self.__Ne + 1)]
        for i in range(self.__Ne - 1):
            en[i + 1] = self.__energy[i][0] + (self.__energy[i + 1][0] - self.__energy[i][0]) / 2
        # Emax in last bin ---!
        en[self.__Ne] = self.__energy[self.__Ne - 1][0] + (self.__energy[self.__Ne - 1][0] - en[self.__Ne - 1])
        # extract spectrum if required ---!
        if self.extract_spectrum:
            self.__extractSpectrumAndModelXML(source_name=source_name, data_path=data_path, scalefluxfactor=scalefluxfactor)
        if return_bin:
            return tbin_stop
        else:
            return

    # get tbin_stop without extracting template data ---!
    def getTimeBinStop(self):
        '''Gets the last time bin of the template if the observation lasts less than the entire afterglow.'''
        self.__getFitsData()
        # time grid ---!
        t = [0.0 for x in range(self.__Nt + 1)]
        for i in range(self.__Nt - 1):
            t[i + 1] = self.__time[i][0] + (self.__time[i + 1][0] - self.__time[i][0]) / 2
        # tmax in last bin ---!
        t[self.__Nt] = self.__time[self.__Nt - 1][0] + (self.__time[self.__Nt - 1][0] - t[self.__Nt - 1])
        # stop the second after higher tmax ---!
        if self.tmax != None:
            tbin_stop = 1
            for bin in range(len(t)):
                if t[bin] <= self.tmax:
                    tbin_stop += 1
                else:
                    continue
        else:
            raise ValueError('Maximum exposure time (tmax) is larger than the template temporal evolution.')
        return tbin_stop

    # get template bins within GTI ---!
    def getTimeBins(self, GTI, tgrid):
        '''Gets which time bins of the template fall within a give time interval.'''
        tbins = []
        for i in range(len(tgrid)):
            # if tgrid[i] <= GTI[0]+10 and tgrid[i+1] >= GTI[0]-10:
            if tgrid[i] <= GTI[0] and tgrid[i+1] >= GTI[0]:
                tbins.append(i)
                continue
            # if tgrid[i] >= GTI[0]-10 and tgrid[i+1] <= GTI[1]+10:
            elif tgrid[i] >= GTI[0] and tgrid[i+1] <= GTI[1]:
                tbins.append(i)
                continue
            # if tgrid[i] >= GTI[1]-10:
            elif tgrid[i] >= GTI[1]:
                tbins.append(i)
                break
        tbins = sorted(tbins)
        tbins = self.__dropListDuplicates(tbins)
        return tbins

    # ctobssim wrapper ---!
    def run_simulation(self, inobs=None, prefix=None, startindex=None):
        '''Wrapper for ctobssim simulation.'''
        self.input = inobs
        sim = ctools.ctobssim()
        if self.input != None: 
            sim["inobs"] = self.input 
        sim["inmodel"] = self.model
        sim["outevents"] = self.output
        sim["caldb"] = self.caldb
        sim["irf"] = self.irf
        if self.edisp:
            sim["edisp"] = self.edisp
        if prefix != None:
            sim["prefix"] = prefix
        if startindex != None:
            sim["startindex"] = startindex
        sim["ra"] = self.pointing[0]
        sim["dec"] = self.pointing[1]
        sim["rad"] = self.fov
        sim["tmin"] = self.t[0]
        sim["tmax"] = self.t[1]
        sim["emin"] = self.e[0]
        sim["emax"] = self.e[1]
        sim["seed"] = self.seed
        sim["nthreads"] = self.nthreads
        sim["logfile"] = self.output.replace('.fits', '.log')
        sim["debug"] = self.set_debug
        if self.set_log:
            sim.logFileOpen()
        sim.execute()
        return

    # dopr duplicates in list ---!
    def __dropListDuplicates(self, list):
        '''Drops duplicate events in list.'''
        new_list = []
        for l in list:
            if l not in new_list:
                new_list.append(l)
        return new_list

    # keep only events within given GTI ---!
    def __dropExceedingEvents(self, hdul, GTI):
        '''Drops events exceeding GTI.'''
        slice_list = []
        times = hdul[1].data.field('TIME')
        for i, t in enumerate(times):
            if t >= GTI[0] and t <= GTI[1]:
                slice_list.append(i)
        return slice_list

    # change from GTI of run to min and max of time events ---!
    def __newGoodTimeIntervals(self, hdul, GTI):
        '''Replaces GTI with min and max time of events.'''
        GTI_new = []
        GTI_new.append(min(hdul[1].data.field('TIME'), key=lambda x: abs(x - GTI[0])))
        GTI_new.append(min(hdul[1].data.field('TIME'), key=lambda x: abs(x - GTI[1])))
        hdul[2].data[0][0] = GTI_new[0]
        hdul[2].data[0][1] = GTI_new[1]
        hdul.flush()
        return

    # reindex rows after sorting ---!
    def __reindexEvents(self, hdul):
        '''Reindexes events.'''
        indexes = hdul[1].data.field(0)
        for i, ind in enumerate(indexes):
            hdul[1].data.field(0)[i] = i + 1
        hdul.flush()
        return

    # sort simulated events by time (TIME) instead of source (MC_ID) ---!
    def __sortEventsByTime(self, hdul, hdr):
        '''Sorts events by time.'''
        data = Table(hdul[1].data)
        data.sort('TIME')
        hdul[1] = fits.BinTableHDU(name='EVENTS', data=data, header=hdr)
        hdul.flush()
        return

    # check GTI and raise error if bad values are passed ---!
    def __checkGTI(self, hdul):
        '''Checks that all events fall within the GTI.'''
        GTI = hdul[2].data[0]
        trange = hdul[1].data.field('TIME')
        if GTI[0] > trange.min() or GTI[1] < trange.max():
            raise ValueError ('Bad GTI values passed to photon list append.')
        return

    # create single photon list from obs list ---!
    def __singlePhotonList(self, sample, filename, GTI, new_GTI=True):
        '''Merge segmented simulations into a single photon list, updating all required header keywords.'''
        sample = sorted(sample)
        n = 0
        for i, f in enumerate(sample):
            with fits.open(f) as hdul:
                if len(hdul[1].data) == 0:
                    continue
                if n == 0:
                    # load header and table ---!
                    hdr1 = hdul[1].header
                    hdr2 = hdul[2].header
                    ext1 = Table(hdul[1].data)
                    ext2 = hdul[2].data
                    n += 1
                else:
                    # update header and append table ---!
                    hdr1['LIVETIME'] += hdul[1].header['LIVETIME']
                    hdr1['ONTIME'] += hdul[1].header['ONTIME']
                    hdr1['TELAPSE'] += hdul[1].header['TELAPSE']
                    hdr1['TSTOP'] = hdul[1].header['TSTOP']
                    hdr1['DATE-END'] = hdul[1].header['DATE-END']
                    hdr1['TIME-END'] = hdul[1].header['TIME-END']
                    ext1 = vstack([ext1, Table(hdul[1].data)])
                hdul.close()
        # create output FITS file empty ---!
        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        hdul.writeto(filename, overwrite=True) if os.path.isfile(filename) else hdul.writeto(filename)
        hdul.close()
        # update FITS file ---!
        with fits.open(filename, mode='update') as hdul:
            hdu1 = fits.BinTableHDU(name='EVENTS', data=ext1, header=hdr1)
            hdu2 = fits.BinTableHDU(name='GTI', data=ext2, header=hdr2)
            hdul.append(hdu1)
            hdul.append(hdu2)
            hdul.flush()
            # sort table by time ---!
            self.__sortEventsByTime(hdul=hdul, hdr=hdr1)
        # manipulate fits ---!
        with fits.open(filename, mode='update') as hdul:
            # drop events exceeding GTI ---!
            slice = self.__dropExceedingEvents(hdul=hdul, GTI=GTI)
            if len(slice) > 0:
                hdul[1].data = hdul[1].data[slice]
            hdul.flush()
            # modify indexes  ---!
            self.__reindexEvents(hdul=hdul)
            # modify GTI ---!
            if new_GTI:
                self.__newGoodTimeIntervals(hdul=hdul, GTI=GTI)
            else:
                hdul[2].data[0][0] = GTI[0]
                hdul[2].data[0][1] = GTI[1]
            hdul.flush()
        return

    # created one FITS table containing all events and GTIs ---!
    def appendEventsSinglePhList(self, GTI=None, new_GTI=False):
        '''From a list of simulations generates a single photon list.'''
        if GTI == None:
            GTI = []
            with fits.open(self.input[0]) as hdul:
                GTI.append(hdul[2].data[0][0])
            with fits.open(self.input[-1]) as hdul:
                GTI.append(hdul[2].data[0][1])
        self.__singlePhotonList(sample=self.input, filename=self.output, GTI=GTI, new_GTI=new_GTI)
        return

    # shift times in template simulation to append background before burst ---!
    def shiftTemplateTime(self, phlist, time_shift):
        '''Shifts events in time.'''
        raise Warning('This method is being fixed.')
        if phlist is str():
            phlist = list(phlist)
        for f in phlist:
            print(f)
            with fits.open(f, mode='update') as hdul:
                # update header ---!
                hdul[1].header['TSTART'] += time_shift
                hdul[1].header['TSTOP'] += time_shift 
                # handle date format to add seconds ---!
                #hdul[1].header['DATE-OBS'] += 
                #hdul[1].header['TIME-OBS'] +=  
                #hdul[1].header['DATE-END'] +=  
                #hdul[1].header['TIME-END'] += 
                # update GTI ---!
                hdul[2].data[0][0] += time_shift
                hdul[2].data[0][1] += time_shift
                # shift events ---!
                if len(hdul[1].data) > 0:
                    times = hdul[1].data.field('TIME')
                    for i, t, in enumerate(times):
                        hdul[1].data.field('TIME')[i] = t + time_shift
                hdul.flush()
        return

    # created a number of observation runs containing all events and GTIs ---!
    def appendEventsMultiPhList(self, max_length=None, last=None, r=True, new_GTI=False):
        '''This method will be deprecated.'''
        exit('This method is outdated')
        n = 1
        sample = []
        singlefile = str(self.output)
        for j in range(int(last / max_length) + 1):
            for i, f in enumerate(self.input):
                with fits.open(f) as hdul:
                    tfirst = hdul[2].data[0][0]
                    tlast = hdul[2].data[0][1]
                    if (tfirst >= max_length * (j) and tlast <= max_length * (j + 1)) or \
                        (tfirst <= max_length * (j) and tlast > max_length * (j)):
                        sample.append(f)
                    elif tlast >= max_length * (j + 1): # or last == max_length * (j + 1):
                        sample.append(f)
                        if n == 1:
                            filename = singlefile.replace('.fits', '_n%03d.fits' % n)
                        else:
                            filename = filename.replace('_n%03d.fits' % (n - 1), '_n%03d.fits' % n)
                        sample = self.__dropListDuplicates(sample)
                        self.__singlePhotonList(sample=sample, filename=filename, GTI=[max_length * (j), max_length * (j + 1)], new_GTI=new_GTI)
                        n += 1
                        drop = len(sample) - 1
                        if drop > 2:
                            self.input = self.input[drop - 1:]
                        sample = [f]
                        break
        if r:
            return n, singlefile
        else:
            return

    def sortObsEvents(self, key='TIME'):
        '''Sorts simulated events by keyword.'''
        with fits.open(self.input, mode='update') as hdul:
            data = Table(hdul[1].data)
            data.sort(key)
            hdr = hdul[1].header
            hdul[1] = fits.BinTableHDU(name='EVENTS', data=data, header=hdr)
            hdul.flush()
            hdul.close()
        with fits.open(self.input, mode='update') as hdul:
            self.__reindexEvents(hdul=hdul)
            hdul.flush()
            hdul.close()
        return