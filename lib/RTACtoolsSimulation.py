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
from astropy import table
from scipy.interpolate import interp1d


class RTACtoolsSimulation():
    '''
    WRITE DOCS
    '''
    def __init__(self):
        # files fields ---!
        self.model, self.template, self.table = (str() for i in range(3))
        self.output, self.input = (str() for i in range(2))
        self.caldb = 'prod2'  # production name in calibration database ---!
        self.irf = 'South_0.5h'  # irf ID name ---!
        # condition control ---!
        self.set_ebl = True  # set/unset EBL absorption feature ---!
        self.extract_spectrum = False  # set/unset spectra extraction feature ---!
        self.plot = False  # option for retrieving plotting values ---!
        self.zfetch = False  # set/unset automatic fetching of redshift ---!
        self.set_debug = False  # set/unset debug mode for ctools ---!
        self.set_log = True  # set/unset logfiles for ctools ---!
        # data ---!
        self.t = [0, 1800]  # time range (s/MJD) ---!
        self.tmax = 1800  # maximum exposure time needed (s) ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.roi = 5  # region of indeterest (deg) ---!
        self.pointing = [83.63, 22.01]  # RA/DEC or GLON/GLAT (deg) ---!
        # ctools miscellaneous ---!
        self.edisp = False  # set/unset edisp
        self.seed = 1  # MC seed ---!
        self.nthreads = 1
        # ebl specifics ---!
        self.z = 0.1  # redshift value ---!
        self.z_ind = 1  # redshift value index ---!
        # fits extension array ---!
        self.__time, self.__energy, self.__spectra, self.__ebl = (float() for i in range(4))

    # open and close the FITS files ---!
    def __openFITS(self):
        hdul = fits.open(self.template)
        return hdul
    def __closeFITS(self, hdul):
        hdul.close()
        return

    # retrive FITS data ---!
    def __getFitsData(self):
        hdul = self.__openFITS()
        self.__energy = np.array(hdul[1].data)
        self.__time = np.array(hdul[2].data)
        self.__spectra = np.array(hdul[3].data)
        if self.set_ebl:
            try:
                self.__ebl = np.array(hdul[4].data)
            except:
                raise IndexError('Template extensions out of range. Unable to load EBL absorbed spectra.')
        self.__closeFITS(hdul)
        return

    # load csv tabl in pandas DataFrame and drop NaN values---!
    def __openCSV(self):
        df = pd.read_csv(self.table)
        df.dropna()
        return df

    # retrive csv data ---!
    def __getEBLfromCSV(self):
        df = self.__openCSV()
        cols = list(df.columns)
        tau_table = np.array(df[cols[self.z_ind]])
        E = np.array(df[cols[0]]) / 1e3  # MeV --> GeV ---!
        return tau_table, E

    # retrive csv temporal bin grid of the template in use and return the necessary slice ---!
    def getTimeSlices(self, GTI, return_bins=False):
        df = self.__openCSV()
        cols = list(df.columns)
        self.__time = np.append(0, np.array(df[cols[1]]))
        bin_start = 0
        bin_stop = 1
        for i in range(len(self.__time)):
            if self.__time[i] < GTI[0]:
                bin_start += 1
                continue
            elif self.__time[i] > GTI[1] or self.__time[i] == GTI[1]:
                self.__time[i] = GTI[1]
                bin_stop += i
                break
        if bin_stop <= self.__Nt:
            time_slice = slice(bin_start, bin_stop + 1)
        else:
            time_slice = slice(bin_start, bin_stop)
        if not time_slice:
            raise ValueError('Invalid GTI: cannot extract time slices')
        if not return_bins:
            return self.__time[time_slice]
        else:
            return self.__time[time_slice], bin_start, bin_stop

    # compute the EBL absorption ---!
    def __addEBL(self, unit='MeV'):
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
        for i in range(len(self.__time)):
            for j in range(len(self.__energy)):
                self.__ebl[i][j] = self.__spectra[i][j] * np.exp(-tau[j])

        if self.plot:
            return E, tau_table, self.__energy, tau
        else:
            return

    # retrive redshift, find nearest column then access its index ---!
    def __zfetch(self):
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
        hdul = self.__openFITS()
        if self.zfetch:
            self.__zfetch()
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
        if self.plot:
            return x, y, x2, y2
        else:
            return

    # extract template spectra, create xml model files and time slices csv file ---!
    def __extractSpectrumAndModelXML(self, source_name, time_slice_name='time_slices.csv', data_path=None):
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
        
            with open(filename, 'a+') as f:
                for j in range(self.__Ne):
                    # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
                    f.write(str(self.__energy[j][0] * 1000.0) + ' ' + str(self.__spectra[i][j] / 1000.0) + "\n")
            # write bin models ---!
            os.system('cp ' + str(self.model) + ' ' + str(os.path.join(data_path, f'{source_name}_tbin{i:02d}.xml')))
            s = open(os.path.join(data_path, f'{source_name}_tbin{i:02d}.xml')).read()
            s = s.replace('data/spec', f'spec_tbin{i:02d}')
            with open(os.path.join(data_path, f'{source_name}_tbin{i:02d}.xml'), 'w') as f:
                f.write(s)
        return

    # read template and return tbin_stop containing necessary exposure time coverage ---!
    def loadTemplate(self, source_name, return_bin=False, data_path=None):
        self.__getFitsData()
        self.__Nt = len(self.__time)
        self.__Ne = len(self.__energy)

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
            raise ValueError('Total exposure time longer than template''s temporal evolution.')

        # energy grid ---!
        en = [1.0 for x in range(self.__Ne + 1)]
        for i in range(self.__Ne - 1):
            en[i + 1] = self.__energy[i][0] + (self.__energy[i + 1][0] - self.__energy[i][0]) / 2
        # Emax in last bin ---!
        en[self.__Ne] = self.__energy[self.__Ne - 1][0] + (self.__energy[self.__Ne - 1][0] - en[self.__Ne - 1])

        # extract spectrum if required ---!
        if self.extract_spectrum:
            self.__extractSpectrumAndModelXML(source_name=source_name, data_path=data_path)
        if return_bin:
            return tbin_stop
        else:
            return

    # get tbin_stop without loading the template ---!
    def getTimeBinStop(self):
        self.__getFitsData()
        self.__Nt = len(self.__time)
        self.__Ne = len(self.__energy)

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
        sim["rad"] = self.roi
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
        new_list = []
        for l in list:
            if l not in new_list:
                new_list.append(l)
        return new_list

    # keep only the GTI rows plus a buffering margin ---!
    def __dropExceedingEvents(self, hdul, GTI):
        slice_list = []
        times = hdul[1].data.field('TIME')
        for i, t in enumerate(times):
            if t >= GTI[0] and t <= GTI[1]:
                slice_list.append(i)
        return slice_list

    # change from GTI of run to min and max of time events ---!
    def __newGoodTimeIntervals(self, hdul, GTI):
        GTI_new = []
        GTI_new.append(min(hdul[1].data.field('TIME'), key=lambda x: abs(x - GTI[0])))
        GTI_new.append(min(hdul[1].data.field('TIME'), key=lambda x: abs(x - GTI[1])))
        hdul[2].data[0][0] = GTI_new[0]
        hdul[2].data[0][1] = GTI_new[1]
        hdul.flush()
        return

    # reindex rows after sorting ---!
    def __reindexEvents(self, hdul):
        indexes = hdul[1].data.field(0)
        for i, ind in enumerate(indexes):
            hdul[1].data.field(0)[i] = i + 1
        hdul.flush()
        return

    # sort simulated events by time (TIME) instead of source (MC_ID) ---!
    def __sortEventsByTime(self, hdul, hdr):
        data = table.Table(hdul[1].data)
        data.sort('TIME')
        hdul[1] = fits.BinTableHDU(name='EVENTS', data=data, header=hdr)
        hdul.flush()
        return

    # create single photon list from obs list ---!
    def __singlePhotonList(self, sample, filename, GTI, new_GTI=False):
        sample = sorted(sample)
        n = 0
        for i, f in enumerate(sample):
            with fits.open(f) as hdul:
                if len(hdul[1].data) == 0:
                    continue
                if n == 0:
                    h1 = hdul[1].header
                    h2 = hdul[2].header
                    ext1 = hdul[1].data
                    ext2 = hdul[2].data
                    n += 1
                else:
                    ext1 = np.append(ext1, hdul[1].data)
        # create output FITS file empty ---!
        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        hdul.writeto(filename, overwrite=True) if os.path.isfile(filename) else hdul.writeto(filename)
        hdul.close()
        # update FITS file ---!
        with fits.open(filename, mode='update') as hdul:
            hdu1 = fits.BinTableHDU(name='EVENTS', data=ext1, header=h1)
            hdu2 = fits.BinTableHDU(name='GTI', data=ext2, header=h2)
            hdul.append(hdu1)
            hdul.append(hdu2)
            hdul.flush()
            # sort table by time ---!
            self.__sortEventsByTime(hdul=hdul, hdr=h1)
        # manipulate fits ---!
        with fits.open(filename, mode='update') as hdul:
            # drop events exceeding GTI ---!
            # slice = self.__dropExceedingEvents(hdul=hdul, GTI=GTI)
            # if len(slice) > 0:
            #     hdul[1].data = hdul[1].data[slice]
            # hdul.flush()
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
        if GTI == None:
            GTI = []
            with fits.open(self.input[0]) as hdul:
                GTI.append(hdul[2].data[0][0])
            with fits.open(self.input[-1]) as hdul:
                GTI.append(hdul[2].data[0][1])
        self.__singlePhotonList(sample=self.input, filename=self.output, GTI=GTI, new_GTI=new_GTI)
        return

    # create one fits table appending source and bkg within GTI ---!
    def appendBkg(self, phlist, bkg, GTI, new_GTI):
        with fits.open(bkg, mode='update') as hdul:
            # fix GTI ---!
            hdul[2].data[0][0] = GTI[0]
            hdul[2].data[0][1] = GTI[1]
            hdul.flush()
            times = hdul[1].data.field('TIME')
            for i, t, in enumerate(times):
                hdul[1].data.field('TIME')[i] = t + GTI[0]
            hdul.flush()
        self.__singlePhotonList(sample=[phlist, bkg], filename=phlist, GTI=GTI, new_GTI=new_GTI)
        return

    # shift times in template simulation to append background before burst ---!
    def shiftTemplateTime(self, phlist, time_shift):
        if phlist is str():
            phlist = list(phlist)
        for file in phlist:
            with fits.open(file, mode='update') as hdul:
                hdul[2].data[0][0] += time_shift
                hdul[2].data[0][1] += time_shift
                times = hdul[1].data.field('TIME')
                for i, t, in enumerate(times):
                    hdul[1].data.field('TIME')[i] = t + time_shift
                hdul.flush()
        return

    # created a number of FITS table containing all events and GTIs ---!
    def appendEventsMultiPhList(self, max_length=None, last=None, r=True, new_GTI=False):
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

    # reduce flux of template by given factor ---!
    def __reduceSpectrumNorm(self, data_path=None):
        spec_files = []
        if data_path is None:
            raise ValueError('please specify a valid path')
        # r: root, d: directories, f: files ---!
        for r, d, f in os.walk(data_path):
            for file in f:
                if self.set_ebl:
                    if '.out' in file and 'ebl' in file and 'flux' not in file:
                        spec_files.append(os.path.os.path.join(r, file))
                else:
                    if '.out' in file and 'ebl' not in file and 'flux' not in file:
                        spec_files.append(os.path.os.path.join(r, file))

        spec_files.sort()
        # new files with relative suffix ---!
        for i in range(len(spec_files)):
            if self.set_ebl:
                new_file = spec_files[i].replace('spec_ebl_tbin', 'spec_ebl_flux%d_tbin' %self.factor)
            else:
                new_file = spec_files[i].replace('spec_tbin', 'spec_flux%d_tbin' %self.factor)
            if os.path.isfile(new_file):
                os.remove(new_file)
            # modify by a given factor ---!
            with open(spec_files[i], 'r') as input, open(new_file, 'w+') as output:
                df = pd.read_csv(input, sep=' ', header=None)
                df.iloc[:,1] = df.iloc[:,1].apply(lambda x: float(x)/self.factor)
                df.to_csv(path_or_buf=output, sep=' ', index=False, header=None)
        return

    # replace path/to/spectrum/file.out in the xml model file ---!
    def __replaceSpectrumFile(self, data_path=None):
        xml_files = []
        if data_path is None:
            raise ValueError('please specify a valid path')
        # r: root, d: directories, f: files ---!
        for r, d, f in os.walk(data_path):
            for file in f:
                if self.set_ebl:
                    if '.xml' in file and 'ebl' in file and 'flux' not in file:
                        xml_files.append(os.path.os.path.join(r, file))
                else:
                    if '.xml' in file and 'ebl' not in file and 'flux' not in file:
                        xml_files.append(os.path.os.path.join(r, file))

        xml_files.sort()
        # replace ---!
        for i in range(len(xml_files)):
            if self.set_ebl:
                new_file = xml_files[i].replace('ID000126_ebl_tbin', 'ID000126_ebl_flux%d_tbin' %self.factor)
            else:
                new_file = xml_files[i].replace('ID000126_tbin', 'ID000126_flux%d_tbin' %self.factor)
            if os.path.isfile(new_file):
                os.remove(new_file)
            with open(xml_files[i], 'r') as input, open(new_file, 'w+') as output:
                content = input.read()
                if self.set_ebl:
                    content = content.replace('spec_ebl_tbin', 'spec_ebl_flux%d_tbin' %self.factor)
                else:
                    content = content.replace('spec_tbin', 'spec_flux%d_tbin' %self.factor)
                output.write(content)
        return

    # execute the flux reduction and consequent substitution of files ---!
    def makeFainter(self, data_path=None):
        self.__reduceSpectrumNorm(data_path=data_path)
        self.__replaceSpectrumFile(data_path=data_path)
        return
