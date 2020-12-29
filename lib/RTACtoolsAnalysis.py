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
import numpy as np

class RTACtoolsAnalysis() :
    '''
    WRITE DOCS
    '''
    def __init__(self):
        # files fields ---!
        self.model = str() 
        self.output, self.input = (str() for i in range(2))
        self.caldb = 'prod2'  # production name in calibration database ---!
        self.irf = 'South_0.5h'  # irf ID name ---!
        # condition control ---!
        self.debug = False  # set/unset debug mode for ctools ---!
        self.if_log = True  # set/unset logfiles for ctools ---!
        # data fields ---!
        self.t = [0, 1800]  # time range (s/MJD) ---!
        self.tmax = 1800  # maximum exposure time needed (s) ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.roi = 5  # region of indeterest (deg) ---!
        self.pointing = [83.63, 22.01]  # RA/DEC or GLON/GLAT (deg) ---!
        self.target = [83.63, 22.51]  # RA/DEC or GLON/GLAT (deg) ---!
        self.sigma = 5  # Gaussian significance (sigmas) ---!
        self.max_src = 10  # Max number of candidates to list during blind-detection ---!
        # ctools miscellaneous ---!
        self.seed = 1  # MC seed ---!
        self.usepnt = True  # use pointing coordinates
        self.refit = False  # refit after initial fit
        self.stats = 'DEFAULT'  # statistics for likelihood fit
        self.accuracy = 0.005  # max like accuracy
        self.max_iter = 50  # max iteration max like
        self.coord_sys = 'CEL'  # coordinate system <CEL|GAL> ---!
        self.sky_subtraction = 'IRF'  # skymap subtraction type <NONE|IRF|RING> ---!
        self.inexclusion = 'NONE'  # exlude sky regions <NONE/file> ---!
        self.bkg_type = 'irf'  # background model <Irf|Aeff|Racc> ---!
        self.src_type = 'POINT'  # source model type ---!
        self.src_name = 'Src001'  # name of source of interest ---!
        self.exclrad = 0.5  # radius around candidate to exclude from further search ---!
        self.corr_kern = 'GAUSSIAN'  # smoothing type ---!
        self.corr_rad = 0.1  # radius for skymap smoothing ---!
        self.fit_shape = False  # enable shape fitting
        self.fit_post = False  # enable position fitting
        self.sgmrange = [0, 10]  # range of gaussian sigmas ---!
        self.confidence = 0.95  # confidence level (%) ---!
        self.eref = 1  # energy reference for flux computation ---!
        self.sens_type = 'Differential'  # sensitivity type <Integral|Differential> ---!
        self.nthreads = 1
        self.stack = False

    # ctselect wrapper ---!
    def run_selection(self, prefix=None):
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
    def run_skymap(self, wbin=0.02, roi_factor=2/np.sqrt(2)):
        nbin = int(self.roi * roi_factor / wbin)
        skymap = ctools.ctskymap()
        skymap['inobs'] = self.input
        skymap['outmap'] = self.output
        skymap['irf'] = self.irf
        skymap['caldb'] = self.caldb
        skymap['emin'] = self.e[0]
        skymap['emax'] = self.e[1]
        skymap['usepnt'] = self.usepnt
        skymap['nxpix'] = nbin
        skymap['nypix'] = nbin
        skymap['binsz'] = wbin
        skymap['coordsys'] = self.coord_sys.upper()
        skymap['proj'] = 'CAR'
        skymap['bkgsubtract'] = self.sky_subtraction.upper()
        skymap['inexclusion'] = self.inexclusion
        # skymap["nthreads"] = self.nthreads
        skymap['logfile'] = self.output.replace('.fits', '.log')
        skymap['debug'] = self.debug
        if self.if_log:
            skymap.logFileOpen()
        skymap.execute()
        return

    # cssrcdetect wrapper ---!
    def run_blindsearch(self) :
        detection = cscripts.cssrcdetect()
        detection['inmap'] = self.input
        detection['outmodel'] = self.output
        detection['outds9file'] = self.output.replace('xml','reg')
        detection['srcmodel'] = self.src_type.upper()
        detection['bkgmodel'] = self.bkg_type.upper()
        detection['fit_pos'] = self.fit_post
        detection['fit_shape'] = self.fit_shape
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

    # csphagen wrapper ---!
    def run_onoff(self, method='reflected', ebins=10, ebins_alg='LOG', binfile=None, exp=None, use_model_bkg=True):
        onoff = cscripts.csphagen()
        onoff['inobs'] = self.input
        onoff['inmodel'] = self.model
        onoff['outobs'] = self.output
        if '.fits' in self.output:
            onoff['outmodel'] = self.output.replace('.fits', '_model.xml')
        elif '.xml' in self.output:
            onoff['outmodel'] = self.output.replace('.xml', '_model.xml')
        else:
            raise ValueError('output obs must be either .fits or .xml file')
        onoff['caldb'] = self.caldb
        onoff['irf'] = self.irf
        onoff['srcname'] = self.src_name
        onoff['ebinalg'] = ebins_alg 
        onoff['emin'] = self.e[0] 
        onoff['emax'] = self.e[1]
        onoff['enumbins'] = ebins 
        if binfile != None:
            onoff['ebinfile'] = binfile
        if exp != None:
            onoff['ebingamma'] = exp
        onoff['coordsys'] = self.coord_sys 
        if self.coord_sys == 'CEL':
            onoff['ra'] = self.target[0]
            onoff['dec'] = self.target[1]
        else:
            onoff['lon'] = self.target[0] 
            onoff['lat'] = self.target[1]
        onoff['rad'] = self.roi 
        onoff['srcregfile'] = self.output.replace('.xml', '_on.reg')
        onoff['bkgregfile'] = self.output.replace('.xml', '_off.reg')
        onoff['bkgmethod'] = method.upper()
        onoff['use_model_bkg'] = use_model_bkg 
        onoff['maxoffset'] = self.roi - 1
        onoff['stack'] = self.stack 
        onoff['etruemin'] = self.e[0] * 0.2 
        onoff['etruemax'] =  self.e[1] * 1.2
        onoff['etruebins'] = round(ebins * 1.5) 
        onoff["nthreads"] = self.nthreads
        onoff['logfile'] = self.output.replace('.xml', '.log')
        onoff['debug'] = self.debug
        if self.if_log:
            onoff.logFileOpen()
        onoff.execute()
        return

    # ctlike wrapper ---!
    def run_maxlikelihood(self):
        like = ctools.ctlike()
        like['inobs'] = self.input
        like['inmodel'] = self.model
        like['outmodel'] = self.output
        like['caldb'] = self.caldb
        like['irf'] = self.irf
        like['refit'] = self.refit
        like['max_iter'] = self.max_iter
        like['like_accuracy'] = self.accuracy
        like['fix_spat_for_ts'] = True
        like['statistic'] = self.stats
        like["nthreads"] = self.nthreads
        like['logfile'] = self.output.replace('.xml', '.log')
        like['debug'] = self.debug
        if self.if_log:
            like.logFileOpen()
        like.execute()
        return

    # cterror wrapper ---!
    def run_asymerrors(self, asym_errors):
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
    def run_uplim(self):
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
    def run_lightcurve(self, nbins=20, bin_type='LIN'):
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

    # cssens wrapper ---!
    def run_sensitivity(self, bins=20, wbin=0.05, enumbins=0):
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

