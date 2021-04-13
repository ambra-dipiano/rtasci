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
import cscripts
from astropy.io import fits
from astropy.table import Table, vstack
import healpy as hp
import numpy as np
import os.path
import pandas as pd
from scipy.interpolate import interp1d
import untangle
import csv
import re
import subprocess
import xml.etree.ElementTree as ET

# configure paths (absolute) ---!
def xmlConfig(cfgfile='/config.xml') :
  # load configuration file ---!
  cfgfile = os.path.dirname(__file__)+str(cfgfile)
  # create configuration dictionary ---!
  with open(cfgfile) as f:
    cfg = untangle.parse(f.read())
  return cfg.config

# analysis dof ---!
def getDof(cfgfile='/config.xml'):
  cfg = xmlConfig(cfgfile)
  if type(cfg.xml.src.free) is list:
    src = len(cfg.xml.src.free)
  else:
    src = len([cfg.xml.src.free])
  bkg = len(cfg.xml.bkg.free)
  dof = src
  return dof, bkg+src, bkg

# true source coords from FITS ---!
def getTrueCoords(fits_file):
  with fits.open(fits_file) as hdul:
    ra = hdul[0].header['RA']
    dec = hdul[0].header['DEC']
  return (ra, dec)

# telescope pointing, either adding off-axis to true coords or as alert probability map peak coords ---!
def getPointing(fits_file, merge_map=None):
  true_coord = getTrueCoords(fits_file)
  if merge_map==None:
    offaxis = (np.random.normal(0,5,1), np.random.normal(0,5,1))
    pointing = (true_coord[0] + offaxis[0], true_coord[1] + offaxis[1])
  else:
    # open and handle map ---!
    map = hp.read_map(merge_map, dtype=None)
    pixels = len(map)
    axis = hp.npix2nside(pixels)
    # search max prob coords ---!
    pmax = np.argmax(map)
    theta, phi = hp.pix2ang(axis, pmax)
    pointing = (np.rad2deg(phi), np.rad2deg(0.5 * np.pi - theta))
    offaxis = (pointing[0] - true_coord[0], pointing[1] - true_coord[1])
  return true_coord, pointing, offaxis

# retrieve telescope pointing coordinates from alert probability map ---!
def getPointingAlert(merge_map=None):
  # load map ---!
  map = hp.read_map(merge_map)
  pixels = len(map)
  axis = hp.npix2nside(pixels)
  # search max prob coords ---!
  pmax = np.argmax(map)
  theta, phi = hp.pix2ang(axis, pmax)
  pointing = (np.rad2deg(phi), np.rad2deg(0.5 * np.pi - theta))
  return pointing

# checks if a trial ID is already existing within a data file ---!
def checkTrialId(file, id):
  with open(file=file) as f:
    df = pd.read_csv(f)
    cols = list(df.columns)
    ids = df[cols[0]]
  if id in list(ids):
    skip = True
  else:
    skip= False
  return skip

# --------------------------------- CLASS xml CONFIGURATION --------------------------------- !!!

class ConfigureXml() :
  '''
  This class handles the configuration of absolute paths for the analysis.
  '''
  def __init__(self, cfg) :
    self.__initPath(cfg)

  # initialise paths ---!
  def __initPath(self, cfg) :
    self.__cfg = cfg
    self.__root = self.__cfg.dir.root['path']
    self.__workdir = self.__cfg.dir.workdir['path']
    self.__runpath = self.__cfg.dir.runpath['path']
    self.__datapath = self.__cfg.dir.datapath['path']
    self.__simpath = self.__cfg.dir.simpath['path']
    self.__selectpath = self.__cfg.dir.selectpath['path']
    self.__detpath = self.__cfg.dir.detpath['path']
    self.__csvpath = self.__cfg.dir.csvpath['path']
    self.__pngpath = self.__cfg.dir.pngpath['path']

  # check the existance of the directory and create if missing ---!
  def __checkDir(self, dir):
    isdir = os.path.isdir(dir)
    return isdir
  def __makeDir(self, dir):
    if not self.__checkDir(dir=dir):
      os.mkdir(dir)

  # get root dir ---!
  def getRootDir(self):
    return self.__root

  # working directory ---!
  def getWorkingDir(self):
    self.__makeDir(self.__workdir.replace('${root}', self.__root))
    return self.__workdir.replace('${root}', self.__root)
  def setWorkingDir(self, workingDir):
    self.__workdir = workingDir
    self.__makeDir(workingDir)

  # directory containing runs ---!
  def getRunDir(self):
    self.__makeDir(self.__runpath.replace('${workdir}', self.getWorkingDir()))
    return self.__runpath.replace('${workdir}', self.getWorkingDir())
  def setRunDir(self, runDir):
    self.__runpath = runDir
    self.__makeDir(runDir)

  # directory storing template data ---!
  def getDataDir(self):
    self.__makeDir(self.__datapath.replace('${runpath}', self.getRunDir()))
    return self.__datapath.replace('${runpath}', self.getRunDir())
  def setDataDir(self, dataDir):
    self.__runpath = dataDir
    self.__makeDir(dataDir)

  # target directory for simulations ---!
  def getSimDir(self):
    self.__makeDir(self.__simpath.replace('${runpath}', self.getRunDir()))
    return self.__simpath.replace('${runpath}', self.getRunDir())
  def setSimDir(self, simDir):
    self.__runpath = simDir
    self.__makeDir(simDir)

  # target directory for selections ---!
  def getSelectDir(self):
    self.__makeDir(self.__selectpath.replace('${runpath}', self.getRunDir()))
    return self.__selectpath.replace('${runpath}', self.getRunDir())
  def setSelectDir(self, selectDir):
    self.__runpath = selectDir
    self.__makeDir(selectDir)

  # target directory for pipeline products ---!
  def getDetDir(self):
    self.__makeDir(self.__detpath.replace('${runpath}', self.getRunDir()))
    return self.__detpath.replace('${runpath}', self.getRunDir())
  def setDetDir(self, detDir):
    self.__runpath = detDir
    self.__makeDir(detDir)

  # target directory for output tables ---!
  def getCsvDir(self):
    self.__makeDir(self.__csvpath.replace('${runpath}', self.getRunDir()))
    return self.__csvpath.replace('${runpath}', self.getRunDir())
  def setCsvDir(self, csvDir):
    self.__runpath = csvDir
    self.__makeDir(csvDir)

  # target directory for images ---!
  def getPngDir(self):
    self.__makeDir(self.__pngpath.replace('${workdir}', self.getWorkingDir()))
    return self.__pngpath.replace('${workdir}', self.getWorkingDir())
  def setPngDir(self, pngDir):
    self.__pngpath = pngDir
    self.__makeDir(pngDir)

# --------------------------------- CLASS ANALYSIS --------------------------------- !!!

class Analysis() :
  '''
  This class contains a mixture of ctools wrapper and pipeline methods. The former are used to easily access and set
  ctools while the latter handles all the analysis necessities: from handling the EBL absorption to degrading the IRFs,
  from extracting spectra to reading the template time bins needed for the simulations. Each public field (self.field)
  can be set accordingly to the user needs from a python script, while private fields (self.__field) is for internal
  usage. Equally, public methods (methodName()) can be invoked within a python script once the class is instanced while
  private methods (__methodName()) should only be used within the class itself.
  '''
  def __init__(self, cfgfile='/config.xml'):
    # location of ctools ---!
    self.__CTOOLS = os.environ.get('CTOOLS')
    # path initial configuration ---!
    self.__cfg = xmlConfig(cfgfile)
    self.__p = ConfigureXml(self.__cfg)
    # files fields ---!
    self.model, self.template, self.table, self.sensCsv = (str() for i in range(4))
    self.output, self.input = (str() for i in range(2))
    self.caldb = 'prod2'  # production name in calibration database ---!
    self.irf = 'South_0.5h'  # irf ID name ---!
    # condition control ---!
    self.if_ebl = True  # set/unset EBL absorption feature ---!
    self.extract_spec = False  # set/unset spectra extraction feature ---!
    self.plot = False  # option for retrieving plotting values ---!
    self.zfetch = False  # set/unset automatic fetching of redshift ---!
    self.debug = False  # set/unset debug mode for ctools ---!
    self.if_log = True  # set/unset logfiles for ctools ---!
    # data fields ---!
    self.t = [0, 1800]  # time range (s/MJD) ---!
    self.tmax = 1800  # maximum exposure time needed (s) ---!
    self.e = [0.03, 150.0]  # energy range (TeV) ---!
    self.roi = 5  # region of indeterest (deg) ---!
    self.pointing = [83.63, 22.01]  # RA/DEC or GLON/GLAT (deg) ---!
    self.sigma = 5  # Gaussian significance (sigmas) ---!
    self.max_src = 10  # Max number of candidates to list during blind-detection ---!
    # ctools miscellaneous ---!
    self.seed = 1  # MC seed ---!
    self.coord_sys = 'CEL'  # coordinate system <CEL|GAL> ---!
    self.sky_subtraction = 'IRF'  # skymap subtraction type <NONE|IRF|RING> ---!
    self.bkg_type = 'irf'  # background model <Irf|Aeff|Racc> ---!
    self.src_type = 'POINT'  # source model type ---!
    self.src_name = 'Src001'  # name of source of interest ---!
    self.exclrad = 0.5  # radius around candidate to exclude from further search ---!
    self.corr_kern = 'GAUSSIAN'  # smoothing type ---!
    self.corr_rad = 0.1  # radius for skymap smoothing ---!
    self.sgmrange = [0, 10]  # range of gaussian sigmas ---!
    self.confidence = 0.95  # confidence level (%) ---!
    self.eref = 1  # energy reference for flux computation ---!
    self.sens_type = 'Differential'  # sensitivity type <Integral|Differential> ---!
    self.nthreads = 1
    # ebl specifics ---!
    self.z = 0.1  # redshift value ---!
    self.z_ind = 1  # redshift value index ---!
    # irf degradation & flux reduction ---!
    self.factor = 2
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
    if self.if_ebl:
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
  def __getCsvData(self):
    df = self.__openCSV()
    cols = list(df.columns)
    tau_gilmore = np.array(df[cols[self.z_ind]])
    E = np.array(df[cols[0]]) / 1e3  # MeV --> GeV ---!
    return tau_gilmore, E

  def getTimeSlices(self, GTI, return_bins=False):
    self.__getFitsData()
    df = self.__openCSV()
    cols = list(df.columns)
    self.__time = np.append(0, np.array(df[cols[1]]))
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
      time_slice = slice(bin_start, bin_stop)
    else:
      time_slice = slice(bin_start, bin_stop-1)
    if not time_slice:
      raise ValueError('Invalid GTI: cannot extract time slices')
    tgrid = self.__time[time_slice]
    tgrid[0] = GTI[0]
    if not return_bins:
      return tgrid
    else:
      return tgrid, bin_start, bin_stop

  # compute the EBL absorption ---!
  def __addEbl(self, unit='MeV'):
    self.__getFitsData()
    tau_gilmore, E = self.__getCsvData()
    if unit == 'GeV':
      E *= 1e3
    elif unit == 'TeV':
      E *= 1e6
    # interpolate linearly handling NaNs/inf/zeroes ---!
    with np.errstate(invalid='raise'):
      interp = interp1d(E, tau_gilmore, bounds_error=False)
    tau = np.array(interp(self.__energy))
    self.__ebl = np.empty_like(self.__spectra)
    # compute absorption ---!
    for i in range(len(self.__time)):
      for j in range(len(self.__energy)):
        self.__ebl[i][j] = self.__spectra[i][j] * np.exp(-tau[j])

    if self.plot:
      return E, tau_gilmore, self.__energy, tau
    else:
      return

  # retrive the redshift and evaluate which table column is the nearest, then access its index ---!
  def __zfetch(self):
    hdul = self.__openFITS()
    # fetch z from the template and chose the table column with min distance from it ---!
    z = hdul[0].header['REDSHIFT']
    with open(self.table, 'r') as f:
      reader = csv.reader(f)
      hdr = next(reader)
    zlist = []
    # load only the redshift columns ---!
    for el in hdr:
      zlist.append(re.sub('[^0-9,.]', '', el))
    zlist.remove('')
    zlist = [float(i) for i in zlist]
    # find nearest ---!
    self.z = min(zlist, key=lambda x:abs(x-z))
    self.z_ind = zlist.index(self.z) +1
    return

  # add EBL extension to a FITS template ---!
  def fitsEbl(self, template_ebl, ext_name='EBL Gilmore', unit='MeV'):
    hdul = self.__openFITS()
    if self.zfetch:
      self.__zfetch()
    if self.plot:
      x, y, x2, y2 = self.__addEbl(unit='MeV')
    else:
      self.__addEbl(unit='MeV')
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
  def __extractSpec(self):
    # time slices table ---!
    table = self.__p.getDataDir() + 'time_slices.csv'
    if os.path.isfile(table):
      os.remove(table)
    with open(table, 'w+') as tab:
      tab.write('#bin,tmax_bin')

    # spectra and models ---!
    for i in range(self.__Nt):
      if self.if_ebl:
        filename = self.__p.getDataDir() + 'spec_ebl_tbin%02d.out' % i
      else:
        filename = self.__p.getDataDir() + 'spec_tbin%02d.out' % i
      if os.path.isfile(filename):
        os.remove(filename)

      # time slices table ---!
      with open(table, 'a') as tab:
        tab.write('\n' + str(i) + ', ' + str(self.__time[i][0]))

      # ebl ---!
      if self.if_ebl:
        with open(filename, 'a+') as f:
          for j in range(self.__Ne):
            # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
            if self.__ebl is not None:
              f.write(str(self.__energy[j][0] * 1000) + ' ' + str(self.__ebl[i][j] / 1000) + "\n")
        # write bin models ---!
        os.system('cp ' + str(self.model) + ' ' + str(self.__p.getDataDir()) + 'run0406_ID000126_ebl_tbin%02d.xml' % i)
        s = open(self.__p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i).read()
        s = s.replace('data/spec', 'spec_ebl_tbin%02d' % i)
        with open(self.__p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i, 'w') as f:
          f.write(s)
      # no ebl ---!
      else:
        with open(filename, 'a+') as f:
          for j in range(self.__Ne):
            # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
            f.write(str(self.__energy[j][0] * 1000.0) + ' ' + str(self.__spectra[i][j] / 1000.0) + "\n")
        # write bin models ---!
        os.system('cp ' + str(self.model) + ' ' + str(self.__p.getDataDir()) + 'run0406_ID000126_tbin%02d.xml' % i)
        s = open(self.__p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i).read()
        s = s.replace('data/spec', 'spec_tbin%02d' % i)
        with open(self.__p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i, 'w') as f:
          f.write(s)
    return

  # read template and return tbin_stop containing necessary exposure time coverage ---!
  def loadTemplate(self) :
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
    if self.tmax != None :
      tbin_stop = 1
      for bin in range(len(t)) :
        if t[bin] <= self.tmax :
          tbin_stop += 1
        else :
          continue
    else :
      raise ValueError('Maximum exposure time (tmax) is larger than the template temporal evolution.')

    # energy grid ---!
    en = [1.0 for x in range(self.__Ne + 1)]
    for i in range(self.__Ne - 1):
      en[i + 1] = self.__energy[i][0] + (self.__energy[i + 1][0] - self.__energy[i][0]) / 2
    # Emax in last bin ---!
    en[self.__Ne] = self.__energy[self.__Ne - 1][0] + (self.__energy[self.__Ne - 1][0] - en[self.__Ne - 1])

    # extract spectrum if required ---!
    if self.extract_spec:
      self.__extractSpec()
    return tbin_stop

  # get tbin_stop without loading the template ---!
  def getTimeBinStop(self) :
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
    if self.tmax != None :
      tbin_stop = 1
      for bin in range(len(t)) :
        if t[bin] <= self.tmax :
          tbin_stop += 1
        else :
          continue
    else :
      raise ValueError('Maximum exposure time (tmax) is larger than the template temporal evolution.')
    return tbin_stop

  # get template bins within GTI ---!
  def getTimeBins(self, GTI, tgrid):
    tbins = []
    for i in range(len(tgrid)):
      if tgrid[i] <= GTI[0]+10 and tgrid[i+1] >= GTI[0]-10:
      # if tgrid[i] <= GTI[0] and tgrid[i+1] >= GTI[0]:
        tbins.append(i)
        continue
      if tgrid[i] >= GTI[0]-10 and tgrid[i+1] <= GTI[1]+10:
      # if tgrid[i] >= GTI[0] and tgrid[i+1] <= GTI[1]:
        tbins.append(i)
        continue
      if tgrid[i] >= GTI[1]-10:
      # if tgrid[i] >= GTI[1]:
        tbins.append(i)
        break

    tbins = sorted(tbins)
    tbins = self.__dropListDuplicates(tbins)
    return tbins

  # ctobssim wrapper ---!
  def eventSim(self) :
    sim = ctools.ctobssim()
    sim["inmodel"] = self.model
    sim["outevents"] = self.output
    sim["caldb"] = self.caldb
    sim["irf"] = self.irf
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
    sim["debug"] = self.debug
    if self.if_log:
      sim.logFileOpen()
    sim.execute()
    return

  # create observation list with gammalib ---!
  def obsList(self, obsname):
    xml = gammalib.GXml()
    obslist = xml.append('observation_list title="observation library"')
    for i in range(len(self.input)):
      obs = obslist.append('observation name="%s" id="%02d" instrument="CTA"' % (obsname, i))
      obs.append('parameter name="EventList" file="%s"' % self.input[i])
    xml.save(self.output)
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
    data = Table(hdul[1].data)
    data.sort('TIME')
    hdul[1] = fits.BinTableHDU(name='EVENTS', data=data, header=hdr)
    hdul.flush()
    return

  # check GTI and raise error if bad values are passed
  def __checkGTI(self, hdul):
    GTI = hdul[2].data[0]
    trange = hdul[1].data.field('TIME')
    if GTI[0] > trange.min() or GTI[1] < trange.max():
      raise ValueError ('Bad GTI values passed to photon list append.')
    return

  # create single photon list from obs list ---!
  def __singlePhotonList(self, sample, filename, GTI, new_GTI=True):
    sample = sorted(sample)
    n = 0
    for i, f in enumerate(sample):
      with fits.open(f) as hdul:
        if len(hdul[1].data) == 0:
          continue
        if n == 0:
          hdr1 = hdul[1].header
          hdr2 = hdul[2].header
          ext1 = Table(hdul[1].data)
          ext2 = hdul[2].data
          n += 1
        else:
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
    if GTI == None:
      GTI = []
      with fits.open(self.input[0]) as hdul:
        GTI.append(hdul[2].data[0][0])
      with fits.open(self.input[-1]) as hdul:
        GTI.append(hdul[2].data[0][1])
    self.__singlePhotonList(sample=self.input, filename=self.output, GTI=GTI, new_GTI=new_GTI)
    return

  # created a number of FITS table containing all events and GTIs ---!
  def appendEventsMultiPhList(self, max_length=None, last=None, r=True):
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
            self.__singlePhotonList(sample=sample, filename=filename, GTI=[max_length * (j), max_length * (j + 1)])
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

  # ctselect wrapper ---!
  def eventSelect(self, prefix=None):
    selection = ctools.ctselect()
    selection['inobs'] = self.input
    selection['outobs'] = self.output
    selection['usepnt'] = True
    if prefix != None:
      selection['prefix'] = prefix
    selection['rad'] = self.roi
    selection['tmin'] = self.t[0]
    selection['tmax'] = self.t[1]
    selection['emin'] = self.e[0]
    selection['emax'] = self.e[1]
    # selection["nthreads"] = self.nthreads
    selection['logfile'] = self.output.replace('.xml', '.log')
    selection['debug'] = self.debug
    if self.if_log:
      selection.logFileOpen()
    selection.execute()
    return

  # ctskymap wrapper ---!
  def eventSkymap(self, wbin=0.02):
    nbin=int(self.roi*2/np.sqrt(2)/wbin)
    skymap = ctools.ctskymap()
    skymap['inobs'] = self.input
    skymap['outmap'] = self.output
    skymap['irf'] = self.irf
    skymap['caldb'] = self.caldb
    skymap['emin'] = self.e[0]
    skymap['emax'] = self.e[1]
    skymap['usepnt'] = True
    skymap['nxpix'] = nbin
    skymap['nypix'] = nbin
    skymap['binsz'] = wbin
    skymap['coordsys'] = self.coord_sys.upper()
    skymap['proj'] = 'CAR'
    skymap['bkgsubtract'] = self.sky_subtraction.upper()
    # skymap["nthreads"] = self.nthreads
    skymap['logfile'] = self.output.replace('.fits', '.log')
    skymap['debug'] = self.debug
    if self.if_log:
      skymap.logFileOpen()
    skymap.execute()
    return

  # cssrcdetect wrapper ---!
  def runDetection(self) :
    self.detectionXml = '%s' % self.input.replace('_skymap.fits', '_det%ssgm.xml' % self.sigma)
    self.detectionReg = '%s' % self.input.replace('_skymap.fits', '_det%ssgm.reg' % self.sigma)

    detection = cscripts.cssrcdetect()
    detection['inmap'] = self.input
    detection['outmodel'] = self.output
    detection['outds9file'] = self.detectionReg
    detection['srcmodel'] = self.src_type.upper()
    detection['bkgmodel'] = self.bkg_type.upper()
    detection['threshold'] = int(self.sigma)
    detection['maxsrcs'] = self.max_src
    detection['exclrad'] = self.exclrad
    detection['corr_rad'] = self.corr_rad
    detection['corr_kern'] = self.corr_kern.upper()
    detection['logfile'] = self.output.replace('.xml', '.log')
    detection['debug'] = self.debug
    # detection["nthreads"] = self.nthreads
    if self.if_log:
      detection.logFileOpen()
    detection.execute()
    return

  # ctlike wrapper ---!
  def maxLikelihood(self):
    like = ctools.ctlike()
    like['inobs'] = self.input
    like['inmodel'] = self.model
    like['outmodel'] = self.output
    like['caldb'] = self.caldb
    like['irf'] = self.irf
    like['refit'] = True
    like['max_iter'] = 500
    like['fix_spat_for_ts'] = False
    like["nthreads"] = self.nthreads
    like['logfile'] = self.output.replace('.xml', '.log')
    like['debug'] = self.debug
    if self.if_log:
      like.logFileOpen()
    like.execute()
    return

  # cterror wrapper ---!
  def confLevels(self, asym_errors):
    self.confidence_level=[0.6827, 0.9545, 0.9973]
    self.output = []
    for i in range(len(self.confidence_level)):
      self.output.append(asym_errors.replace('_errors', '_%2derr' % (self.confidence_level[i] * 100)))
      if not os.path.isfile(self.output[i]):
        err = ctools.cterror()
        err['inobs'] = self.input
        err['inmodel'] = self.model
        err['srcname'] = self.src_name
        err['outmodel'] = self.output[i]
        err['caldb'] = self.caldb
        err['irf'] = self.irf
        err['confidence'] = self.confidence_level[i]
        err["nthreads"] = self.nthreads
        err['logfile'] = self.output[i].replace('.xml', '.log')
        err['debug'] = self.debug
        if self.if_log:
          err.logFileOpen()
        err.execute()
    return self.output

  # ctulimit wrapper ---!
  def integrFlux(self):
    uplim = ctools.ctulimit()
    uplim['inobs'] = self.input
    uplim['inmodel'] = self.model
    uplim['srcname'] = self.src_name
    uplim['caldb'] = self.caldb
    uplim['irf'] = self.irf
    uplim['confidence'] = self.confidence
    uplim['sigma_min'] = self.sgmrange[0]
    uplim['sigma_max'] = self.sgmrange[1]
    uplim['eref'] = self.eref  # default reference energy for differential limit (in TeV)
    uplim['emin'] = self.e[0]  # default minimum energy for integral flux limit (in TeV)
    uplim['emax'] = self.e[1]  # default maximum energy for integral flux limit (in TeV)
    uplim["nthreads"] = self.nthreads
    uplim['logfile'] = self.model.replace('results.xml', 'flux.log')
    uplim['debug'] = self.debug
    if self.if_log:
      uplim.logFileOpen()
    uplim.execute()
    return

  # ctulimit wrapper ---!
  def ctoolsLightCurve(self, nbins=20, bin_type='LIN'):
    lc = ctools.cslightcrv()
    lc['inobs'] = self.input
    lc['inmodel'] = self.model
    lc['srcname'] = self.src_name
    lc['caldb'] = self.caldb
    lc['irf'] = self.irf
    lc['outfile'] = self.output
    lc['tbinalg'] = bin_type
    lc['tmin'] = self.t[0]
    lc['tmax'] = self.t[1] # default reference energy for differential limit (in TeV)
    lc['tbins'] = nbins
    lc['method'] = '3D'
    lc['emin'] = self.e[0]  # default minimum energy for integral flux limit (in TeV)
    lc['emax'] = self.e[1]  # default maximum energy for integral flux limit (in TeV)
    lc['enumbins'] = 0
    lc['coordsys'] = self.coord_sys
    lc['proj'] = 'CAR'
    lc['xref'] = self.pointing[0]
    lc['yref'] = self.pointing[1]
    lc["nthreads"] = self.nthreads
    lc['logfile'] = self.output.replace('.xml', '.log')
    lc['debug'] = self.debug
    if self.if_log:
      lc.logFileOpen()
    lc.execute()
    return

  # compute integral photon flux for PL model ---!
  def photonFluxPowerLaw(self, gamma, k0, e0=1, unit='TeV'):
    if unit == 'eV':
      conv = 1e-6
    elif unit == 'keV':
      conv = 1e-3
    elif unit == 'MeV':
      conv = 1
    elif unit == 'GeV':
      conv = 1e3
    else:
      conv = 1e6
    e1 = self.e[0] * conv
    e2 = self.e[1] * conv
    delta = gamma + 1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

  # compute integral energy flux for PL model ---!
  def energyFluxPowerLaw(self, gamma, k0, e0=1, unit='TeV'):
    if unit == 'eV':
      conv = 1.60218e-12
    elif unit == 'keV':
      conv = 1.60218e-9
    elif unit == 'MeV':
      conv = 1.60218e-6
    elif unit == 'GeV':
      conv = 1.60218e-3
    else:
      conv = 1.60218
    e1 = self.e[0] * conv
    e2 = self.e[1] * conv
    k0 *= conv
    e0 *= conv
    delta = gamma+1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

  # initialize paths for caldb degradation: directories and files ---!
  def __initCaldbIrf(self):
    nominal_irf =  self.__CTOOLS + '/share/caldb/data/cta/%s/bcf/%s/irf_file.fits' % (self.caldb, self.irf)
    degraded_irf = nominal_irf.replace('prod', 'degr')
    caldb_degr = self.caldb.replace('prod', 'degr')
    folder = self.__CTOOLS + '/share/caldb/data/cta/'
    nominal_cal =  folder + self.caldb
    degraded_cal = folder + caldb_degr
    return folder, nominal_cal, nominal_irf, degraded_cal, degraded_irf

  # updates the degraded caldb index by replacing all "prod" references with "degr" ---!
  def __updateCaldbIndex(self, index):
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
    if not os.path.isdir(degraded_cal):
      os.mkdir(degraded_cal)
    if not os.path.isfile(degraded_cal+'/caldb.indx'):
      os.system('cp %s/caldb.indx %s/caldb.indx' %(nominal_cal, degraded_cal))
      # update caldb.indx file ---!
      self.__updateCaldbIndex(degraded_cal+'/caldb.indx')
    if not os.path.isdir(degraded_cal+'/bcf'):
      os.mkdir(degraded_cal+'/bcf')
    if not os.path.isdir(degraded_cal+'/bcf/'+self.irf):
      os.mkdir(degraded_cal+'/bcf/'+self.irf)
    if os.path.isfile(degraded_irf):
      os.system('rm %s' %degraded_irf)
    if not os.path.isfile(degraded_irf):
      os.system('cp %s %s' %(nominal_irf, degraded_irf))
    return

  # change permission to 777 and ask for password if user id not in idlist param ---!
  def __openPermission(self, path, idlist=(0,1126)):
    if os.geteuid() in idlist:
      subprocess.run(['chmod', '-R', '777', path], check=True)
    else:
      subprocess.run(['sudo', 'chmod', '-R', '777', path], check=True)
    return

  # change permission to 755 and ask for password if user id not in idlist param ---!
  def __closePermission(self, path, idlist=(0,1126)):
    if os.geteuid() in idlist:
      subprocess.run(['chmod', '-R', '755', path], check=True)
    else:
      subprocess.run(['sudo', 'chmod', '-R', '755', path], check=True)
    return

  # degrade Aff by self.factor (either scalar or array of energy-bins dimension) ---!
  def __degrAeff(self, nominal_irf, degraded_irf, r=False):
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
  def __degrBkg(self, nominal_irf, degraded_irf, aeff=True):
    # degrade Aeff (only if True) and get its returns ---!
    if not aeff:
      tmp = self.factor
      self.factor = 1
    aeff_nom, aeff_deg, theta, e_aeff = self.__degrAeff(nominal_irf=nominal_irf, degraded_irf=degraded_irf, r=True)
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
  def degradeIrf(self, bkg=True, aeff=True):
    # initialize ---!
    folder, nominal_cal, nominal_irf, degraded_cal, degraded_irf = self.__initCaldbIrf()
    # open all folder permission ---!
    self.__openPermission(path=folder)
    # create degr caldb path if not existing ---!
    self.__mockNominalCaldb(nominal_cal=nominal_cal, nominal_irf=nominal_irf,
                            degraded_cal=degraded_cal, degraded_irf=degraded_irf)
    # close all folder permission and open only degraded caldb permission ---!
    self.__closePermission(path=folder)
    self.__openPermission(path=degraded_cal)
    # degradation aeff ---!
    if not bkg:
      self.__degrAeff(nominal_irf=nominal_irf, degraded_irf=degraded_irf)
    # degradation bkg counts ---!
    else:
      self.__degrBkg(nominal_irf=nominal_irf, degraded_irf=degraded_irf, aeff=aeff)
    # close degraded caldb permission ---!
    self.__closePermission(degraded_cal)
    # update caldb ---!
    self.caldb = self.caldb.replace('prod', 'degr')
    return

  # cssens wrapper ---!
  def eventSens(self, bins=1, wbin=0.05, enumbins=0):
    sens = cscripts.cssens()
    nbin = int(self.roi / wbin)
    sens['inobs'] = self.input
    sens['inmodel'] = self.model
    sens['srcname'] = self.src_name
    sens['caldb'] = self.caldb
    sens['irf'] = self.irf
    sens['outfile'] = self.output
    sens['duration'] = self.t[1] - self.t[0]
    sens['rad'] = self.roi
    sens['emin'] = self.e[0]
    sens['emax'] = self.e[1]
    sens['bins'] = bins
    if enumbins != 0:
      sens['enumbins'] = enumbins
      sens['npix'] = nbin
      sens['binsz'] = wbin
    sens['sigma'] = self.sigma
    sens['type'] = self.sens_type.capitalize()
    sens["nthreads"] = self.nthreads
    sens['logfile'] = self.output.replace('.csv', '.log')
    sens['debug'] = self.debug
    if self.if_log:
      sens.logFileOpen()
    sens.execute()
    return

  # reduce flux of template by given factor ---!
  def __reduceFluxSpec(self):
    spec_files = []
    # r: root, d: directories, f: files ---!
    for r, d, f in os.walk(self.__p.getDataDir()):
      for file in f:
        if self.if_ebl:
          if '.out' in file and 'ebl' in file and 'flux' not in file:
            spec_files.append(os.path.join(r, file))
        else:
          if '.out' in file and 'ebl' not in file and 'flux' not in file:
            spec_files.append(os.path.join(r, file))

    spec_files.sort()
    # new files with relative suffix ---!
    for i in range(len(spec_files)):
      if self.if_ebl:
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
  def __replaceSpecFile(self):
    xml_files = []
    # r: root, d: directories, f: files ---!
    for r, d, f in os.walk(self.__p.getDataDir()):
      for file in f:
        if self.if_ebl:
          if '.xml' in file and 'ebl' in file and 'flux' not in file:
            xml_files.append(os.path.join(r, file))
        else:
          if '.xml' in file and 'ebl' not in file and 'flux' not in file:
            xml_files.append(os.path.join(r, file))

    xml_files.sort()
    # replace ---!
    for i in range(len(xml_files)):
      if self.if_ebl:
        new_file = xml_files[i].replace('ID000126_ebl_tbin', 'ID000126_ebl_flux%d_tbin' %self.factor)
      else:
        new_file = xml_files[i].replace('ID000126_tbin', 'ID000126_flux%d_tbin' %self.factor)
      if os.path.isfile(new_file):
        os.remove(new_file)
      with open(xml_files[i], 'r') as input, open(new_file, 'w+') as output:
        content = input.read()
        if self.if_ebl:
          content = content.replace('spec_ebl_tbin', 'spec_ebl_flux%d_tbin' %self.factor)
        else:
          content = content.replace('spec_tbin', 'spec_flux%d_tbin' %self.factor)
        output.write(content)
    return

  # execute the flux reduction and consequent substitution of files ---!
  def makeFainter(self):
    self.__reduceFluxSpec()
    self.__replaceSpecFile()
    return

  # returns a random total delay time (slew time + gw latency) within given ranges ---!
  def totalDelay(self, slew=(0,50), gw_latency=(0,36000)):
    tslew = np.random.uniform(slew[0], slew[1], 1)
    tgw = np.random.uniform(gw_latency[0], gw_latency[1])
    delay = tslew + tgw
    return delay

# --------------------------------- CLASS xml HANDLING --------------------------------- !!!

class ManageXml():
  '''
  This class contains all the methods which read, generate and modify xml files, as needed for the analysis.
  They are not comprehensive of all the parameters one could want to access though more could be added according
  to necessity. In the future this class will also handle the pipeline configuration files. Each public field (self.field)
  can be set accordingly to the user needs from a python script, while private fields (self.__field) is for internal
  usage. Equally, public methods (methodName()) can be invoked within a python script once the class is instanced while
  private methods (__methodName()) should only be used within the class itself.
  '''
  def __init__(self, xml, cfgfile='/config.xml'):
    self.__xml = xml
    self.__cfg = xmlConfig(cfgfile)
    self.__p = ConfigureXml(self.__cfg)
    self.file = open(self.__xml)
    self.src_lib = ET.parse(self.file)
    self.root = self.src_lib.getroot()
    self.tsv_list = []
    self.pos = []
    self.err = []
    self.spectral = []
    self.sigma = 5
    self.default_model = True
    self.instr = 'CTA'
    self.bkg_type = 'Irf'
    self.src_att = []
    self.bkg_att = []
    self.tscalc = True
    self.if_cut = False

  # get source element ---!
  def __getSrcObj(self):
    src = self.root.findall('source')
    return src

  # skip node listed in skip element or filters ---!
  def __skipNode(self, cfg):
    src = self.__getSrcObj()
    if src.attrib[cfg.get('idAttribute')] in cfg.get('skip'):
      return True
    for filter in cfg.get('filters'):
      if src.attrib[filter.get('attribute')] == filter.get('value'):
        return True
    if len(cfg.get('selectors')) == 0:
      return False
    for select in cfg.get('selectors'):
      if src.attrib[select.get('attribute')] == select.get('value'):
        return False
    return True

  # get TS values ---!
  def loadTs(self, highest=None):
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          tsv = src.attrib['ts']
          self.tsv_list.append(tsv)
      else:
        if src.attrib['name'] == highest:
          tsv = src.attrib['ts']
          self.tsv_list.append(tsv)
    return self.tsv_list

  # get RA/DEC values ---!
  def loadRaDec(self, highest=None):
    ra_list, dec_list = ([] for i in range(2))
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          ra_list.append(ra)
          dec_list.append(dec)
      else:
        if src.attrib['name'] == highest:
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          ra_list.append(ra)
          dec_list.append(dec)
    self.pos = [ra_list, dec_list]
    return self.pos

  # get Gaussian sigma values ---!
  def loadConfInt(self, highest=None):
    ra_list, dec_list = ([] for i in range(2))
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          ra_list.append(ra)
          dec_list.append(dec)
      else:
        if src.attrib['name'] == highest:
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          ra_list.append(ra)
          dec_list.append(dec)
    self.err = [ra_list, dec_list]
    return self.err

  # get spectral parameter values ---1
  def loadSpectral(self, highest=None):
    index_list, pref_list, pivot_list = ([] for i in range(3))
    if self.if_cut is True :
      cutoff_list = []

    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        if highest == None:
          index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
              src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
          pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
              src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
          pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
              src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
          index_list.append(index)
          pref_list.append(pref)
          pivot_list.append(pivot)
          if self.if_cut is True :
            cutoff = float(src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['scale'])
            cutoff_list.append(cutoff)
        else:
          if src.attrib['name'] == highest:
            index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
            pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
            pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
            index_list.append(index)
            pref_list.append(pref)
            pivot_list.append(pivot)
            if self.if_cut is True:
              cutoff = float(src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['value']) * float(
                  src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['scale'])
              cutoff_list.append(cutoff)

    if self.if_cut is False:
      self.spectral = [index_list, pref_list, pivot_list]
    else:
      self.spectral = [index_list, pref_list, pivot_list, cutoff_list]
    return self.spectral

  # get prefact error values ---!
  def loadPrefError(self, highest=None):
    err_list = []
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          err = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['error']) * float(
              src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
          err_list.append(err)
      else:
        if src.attrib['name'] == highest:
          err = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['error']) * float(
              src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
          err_list.append(err)
    self.err = err_list
    return self.err

  # save xml files ---!
  def __saveXml(self):
    self.src_lib.write(self.__xml, encoding="UTF-8", xml_declaration=True)
    return

  # set a default spectral model and bkg ---!
  def __setModel(self):
    if self.default_model is True:
      att_prefactor = {'name': 'Prefactor', 'scale': '1e-16', 'value': '5.7', 'min': '1e-07', 'max': '1e7', 'free': '1'}
      att_index = {'name': 'Index', 'scale': '-1', 'value': '2.48', 'min': '0', 'max': '5.0', 'free': '1'}
      att_pivot = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1.0', 'min': '1e-07', 'max': '1000.0', 'free': '0'}
      bkg_prefactor = {'name': 'Prefactor', 'scale': '1', 'value': '1.0', 'min': '1e-03', 'max': '1e+3.0', 'free': '1'}
      bkg_index = {'name': 'Index', 'scale': '1.0', 'value': '0.0', 'min': '-5', 'max': '+5.0', 'free': '1'}
      bkg_pivot = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '0'}

      self.src_att = [att_prefactor, att_index, att_pivot]
      self.bkg_att = [bkg_prefactor, bkg_index, bkg_pivot]
      if self.if_cut is True:
        att_cutoff = {'name': 'CutoffEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '1'}
        self.src_att.append(att_cutoff)

      return self.src_att, self.bkg_att
    else:
      pass

  # set tscalc to 1 ---!
  def setTsTrue(self):
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        src.set('tscalc', '1')
    self.__saveXml()
    return

  # modeify the spectral component of candidate list ---!
  def modXml(self, overwrite=True):
    self.__setModel()
    # source ---!
    i = 0
    for src in self.root.findall('source'):
      i += 1
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        src.set('tscalc', '1') if self.tscalc is True else None
        # remove spectral component ---!
        rm = src.find('spectrum')
        src.remove(rm)
        # new spectrum ---!
        if self.if_cut:
          spc = ET.SubElement(src, 'spectrum', attrib={'type': 'ExponentialCutoffPowerLaw'})
        else:
          spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
        spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
        spc.tail = '\n\t\t'.replace('\t', ' ' * 2)
        # new spectral params ---!
        for j in range(len(self.src_att)):
          prm = ET.SubElement(spc, 'parameter', attrib=self.src_att[j])
          if prm.attrib['name'] == 'Prefactor' and i > 1:
            prm.set('value', str(float(prm.attrib['value']) / 2 ** (i - 1)))
          prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < len(self.src_att) else '\n\t\t'.replace('\t', ' ' * 2)
          #spc.insert(j, prm)
      # background ---!
      else:
        # set bkg attributes ---!
        src.set('instrument', '%s' % self.instr.upper()) if self.instr.capitalize() != 'None' else None
        if self.bkg_type.capitalize() == 'Aeff' or self.bkg_type.capitalize() == 'Irf':
          src.set('type', 'CTA%sBackground' % self.bkg_type.capitalize())
        if self.bkg_type.capitalize() == 'Racc':
          src.set('type', 'RadialAcceptance')
        # remove spectral component ---!
        rm = src.find('spectrum')
        src.remove(rm)
        # new bkg spectrum ---!
        spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
        spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
        spc.tail = '\n\t'.replace('\t', ' ' * 2)
        # new bkg params ---!
        for j in range(len(self.bkg_att)):
          prm = ET.SubElement(spc, 'parameter', attrib=self.bkg_att[j])
          prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < len(self.bkg_att) else '\n\t\t'.replace('\t', ' ' * 2)

    # instead of override original xml, save to a new one with suffix "_mod" ---!
    if not overwrite:
      self.__xml = self.__xml.replace('.xml', '_mod.xml')
    self.__saveXml()
    return

  # free and fix parameters for max like computation ---!
  def prmsFreeFix(self):
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        for prm in src.findall('*/parameter'):
          if prm.attrib['name'] not in self.__cfg.xml.bkg.free:
            prm.set('free', '0')
        for free in self.__cfg.xml.src.free:
          src.find('*/parameter[@name="%s"]' % free['prm']).set('free', '1') if free['prm'] != None else None
      else:
        for prm in src.findall('*/parameter'):
          if prm.attrib['name'] not in self.__cfg.xml.bkg.free:
            prm.set('free', '0')
        for free in self.__cfg.xml.bkg.free:
          src.find('*/parameter[@name="%s"]' % free['prm']).set('free', '1') if free['prm'] != None else None

    #self.setTsTrue() if self.tscalc is True else None
    self.__saveXml()
    return

  # sort candidates by their ts value ---!
  def sortSrcTs(self):
    src = self.root.findall("*[@ts]")
    self.root[:-1] = sorted(src, key=lambda el: (el.tag, el.attrib['ts']), reverse=True)
    from_highest = []
    for src in self.root.findall("*[@ts]"):
      from_highest.append(src.attrib['name'])
    self.__saveXml()
    if len(from_highest) == 0:
      from_highest = [None]
    return from_highest

  # close xml ---!
  def closeXml(self):
    self.file.close()
    return

