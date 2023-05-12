# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import ctools
import cscripts
import os
import numpy as np
from astropy.io import fits

# count photometry ctools
def onoff_counts(pha):
    '''This function returns on, off and excess counts given ctools on/off files.'''
    regions = pha.replace('.xml', '_off.reg')
    off = pha.replace('.xml', '_pha_off.fits')
    on = pha.replace('.xml', '_pha_on.fits')
    nreg = 0
    # regions ---!
    with open(regions, 'r') as f:
        for line in f:
            if line.startswith('fk5'):
                nreg += 1
    # on counts ---!
    with fits.open(on) as h:
        oncounts = h[1].data.field('counts')
    onsum = 0
    for val in oncounts:
        onsum += val
    # off counts ---!
    with fits.open(off) as h:
        offcounts = h[1].data.field('counts')
    offsum = 0
    for val in offcounts:
        offsum += val
    offsum = np.sum(offsum)
    # alpha and excess ---!
    alpha = 1.0/nreg
    excess = onsum - alpha * offsum
    return onsum, offsum, excess, alpha

class RTACtoolsAnalysis() :
    '''
    This class contains wrappers for ctools and cscripts tools.
    '''
    def __init__(self, on_ram=False):
        # files fields ---!
        self.__on_ram = on_ram
        self.model = str() 
        self.output, self.input = (str() for i in range(2))
        self.caldb = 'prod2'  # caldb (str) ---!
        self.irf = 'South_0.5h'  # irf (Str) ---!
        # condition control ---!
        self.set_debug = False  # set/unset debug mode for ctools ---!
        self.set_log = True  # set/unset logfiles for ctools ---!
        self.edisp = True  # require enegy dispersion ---!
        # data ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.roi = 5  # region of interest (deg) ---!
        self.t = [0, 1800]  # time range (s/MJD) ---!
        self.pointing = [83.63, 22.01]  # RA/DEC or GLON/GLAT (deg) ---!
        self.target = [83.63, 22.51]  # RA/DEC or GLON/GLAT (deg) ---!
        self.sigma = 3  # Gaussian significance (sigmas) ---!
        self.max_src = 10  # Max number of candidates to list during blind-detection ---!
        # ctools miscellaneous ---!
        self.seed = 1  # MC seed ---!
        self.usepnt = True  # use pointing coordinates
        self.refit = False  # refit after initial fit
        self.stats = 'DEFAULT'  # statistics for likelihood fit
        self.accuracy = 0.005  # max like accuracy
        self.max_iter = 50  # max iteration max like
        self.coord_sys = 'CEL'  # coordinate system <CEL|GAL> ---!
        self.proj = 'CAR'  # projection method ---!
        self.sky_subtraction = 'IRF'  # skymap subtraction type <NONE|IRF|RING> ---!
        self.inexclusion = 'NONE'  # exlude sky regions <NONE/file> ---!
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
        self.stack = False

    # ctselect wrapper ---!
    def run_selection(self, prefix=None):
        '''Wrapper of ctselect.'''
        selection = ctools.ctselect()
        selection['inobs'] = self.input
        selection['outobs'] = self.output
        selection['usepnt'] = self.usepnt
        if prefix != None:
            selection['prefix'] = prefix
        selection['rad'] = self.roi
        selection['tmin'] = self.t[0]
        selection['tmax'] = self.t[1]
        selection['emin'] = self.e[0]
        selection['emax'] = self.e[1]
        # selection["nthreads"] = self.nthreads
        selection['logfile'] = self.output.replace('.xml', '.log')
        selection['debug'] = self.set_debug
        if self.set_log:
            selection.logFileOpen()
        if not self.__on_ram:
            selection.execute()
        else:
            selection.run()
        return

    # ctskymap wrapper ---!
    def run_skymap(self, wbin=0.02, roi_factor=1):
        '''Wrapper of ctskymap.'''
        nbin = int(self.roi*2*roi_factor/wbin)
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
        skymap['debug'] = self.set_debug
        if self.set_log:
            skymap.logFileOpen()
        if not self.__on_ram:
            skymap.execute()
        else:
            skymap.run()
        return

    # cssrcdetect wrapper ---!
    def run_blindsearch(self, fit_pos=False, fit_shape=False):
        '''Wrapper of cssrcdetect.'''
        detection = cscripts.cssrcdetect()
        detection['inmap'] = self.input
        detection['outmodel'] = self.output
        detection['outds9file'] = self.output.replace('xml','reg')
        detection['srcmodel'] = self.src_type.upper()
        detection['bkgmodel'] = self.bkg_type.upper()
        detection['fit_pos'] = fit_pos
        detection['fit_shape'] = fit_shape
        detection['threshold'] = int(self.sigma)
        detection['maxsrcs'] = self.max_src
        detection['exclrad'] = self.exclrad
        detection['corr_rad'] = self.corr_rad
        detection['corr_kern'] = self.corr_kern.upper()
        detection['logfile'] = self.output.replace('.xml', '.log')
        detection['debug'] = self.set_debug
        # detection["nthreads"] = self.nthreads
        if self.set_log:
            detection.logFileOpen()
        if not self.__on_ram:
            detection.execute()
        else:
            detection.run()
        return

    # csphagen wrapper ---!
    def run_onoff(self, method='reflected', prefix='onoff', maxoffset=2.5, radius=0.2, ebins=40, ebins_alg='LOG', binfile=None, exp=None, use_model_bkg=True, etruemin=0.01, etruemax=0.01, etruebins=30, bkgskip=1, bkgmin=2):
        '''Wrapper of csphagen.'''
        onoff = cscripts.csphagen()
        onoff['inobs'] = self.input
        onoff['inmodel'] = self.model
        onoff['prefix'] = prefix
        onoff['outobs'] = self.output
        if '.xml' in self.output:
            onoff['outmodel'] = self.output.replace('.xml', '_model.xml')
        else:
            raise ValueError('output obs must be a .xml file')
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
        onoff['rad'] = radius
        onoff['srcregfile'] = self.output.replace('.xml', '_on.reg')
        onoff['bkgregfile'] = self.output.replace('.xml', '_off.reg')
        onoff['bkgmethod'] = method.upper()
        onoff['use_model_bkg'] = use_model_bkg 
        onoff['maxoffset'] = maxoffset
        onoff['bkgregmin'] = bkgmin
        onoff['bkgregskip'] = bkgskip
        onoff['stack'] = self.stack 
        onoff['etruemin'] = etruemin
        onoff['etruemax'] =  etruemax
        onoff['etruebins'] = etruebins 
        onoff["nthreads"] = self.nthreads
        onoff['logfile'] = self.output.replace('.xml', '.log')
        onoff['debug'] = self.set_debug
        if self.set_log:
            onoff.logFileOpen()
        if not self.__on_ram:
            onoff.execute()
        else:
            onoff.run()
        return

    # ctbin wrapper ---!
    def run_binning(self, prefix='cube_', ebins_alg='LOG', ebins=10, binfile=None, exp=None, nbins=None, wbin=0.02):
        '''Wrapper of ctbin.'''
        if nbins == None:
            nbins = int(self.roi * 2 / wbin)
        bins = ctools.ctbin()
        bins['inobs'] = self.input
        bins['outobs'] = self.output
        bins['stack'] = self.stack
        if prefix != None:
            bins['prefix'] = prefix
        bins['ebinalg'] = ebins_alg
        bins['emin'] = self.e[0]
        bins['emax'] = self.e[1]
        bins['enumbins'] = ebins
        if binfile != None:
            bins['ebinfile'] = binfile 
        if exp != None:
            bins['ebingamma'] = exp
        bins['usepnt'] = self.usepnt
        bins['nxpix'] = nbins
        bins['nypix'] = nbins
        bins['binsz'] = wbin
        bins['coordsys'] = self.coord_sys
        bins['proj'] = self.proj
        if self.usepnt is False:
            bins['xref'] = self.target[0] 
            bins['yref'] = self.target[1] 
        if '.fits' in self.output:
            bins['logfile'] = self.output.replace('.fits', '.log')
        elif '.xml' in self.output:
            bins['logfile'] = self.output.replace('.xml', '.log')
        bins['debug'] = self.set_debug
        if self.set_log:
            bins.logFileOpen()
        if not self.__on_ram:
            bins.execute()
        else:
            bins.run()
        return

    # ctexpcube wrapper ---!
    def run_expcube(self, cube, ebins=10, nbins=None, wbin=0.02, ebin_alg='LOG', ebinfile=None, ebingamma=None, addbounds=False):
        '''Wrapper of ctexpcube.'''
        if nbins == None:
            nbins = int(self.roi * 2 / wbin)
        exp = ctools.ctexpcube()
        exp['inobs'] = self.input
        exp['incube'] = cube
        exp['caldb'] = self.caldb
        exp['irf'] = self.irf
        exp['outcube'] = self.output
        exp['ebinalg'] = ebin_alg
        exp['emin'] = self.e[0]
        exp['emax'] = self.e[1]
        exp['enumbins'] = ebins
        if ebinfile != None:
            exp['ebinfile'] = ebinfile
        if ebingamma != None:
            exp['ebingamma'] = ebingamma
        exp['addbounds'] = addbounds
        exp['usepnt'] = self.usepnt
        exp['nxpix'] = nbins
        exp['nypix'] = nbins
        exp['binsz'] = wbin
        exp['coordsys'] = self.coord_sys
        exp['proj'] = self.proj
        if not self.usepnt:
            exp['xref'] = self.target[0] 
            exp['yref'] = self.target[1] 
        if '.fits' in self.output:
            exp['logfile'] = self.output.replace('.fits', '.log')
        elif '.xml' in self.output:
            exp['logfile'] = self.output.replace('.xml', '.log')
        exp['debug'] = self.set_debug
        if self.set_log:
            exp.logFileOpen()
        if not self.__on_ram:
            exp.execute()
        else: 
            exp.run()
        return
    
    # ctpsfcube wrapper ---!
    def run_psfcube(self, cube, ebins=10, nbins=None, wbin=0.02, amax=0.3, abins=200, ebin_alg='LOG', ebinfile=None, ebingamma=None, addbounds=False):
        '''Wrapper of ctpsfcube.'''
        if nbins == None:
            nbins = int(self.roi * 2 / wbin)
        psf = ctools.ctpsfcube()
        psf['inobs'] = self.input
        psf['incube'] = cube
        psf['caldb'] = self.caldb
        psf['irf'] = self.irf
        psf['outcube'] = self.output
        psf['ebinalg'] = ebin_alg
        psf['emin'] = self.e[0]
        psf['emax'] = self.e[1]
        psf['enumbins'] = ebins
        if ebinfile != None:
            psf['ebinfile'] = ebinfile
        if ebingamma != None:
            psf['ebingamma'] = ebingamma
        psf['addbounds'] = addbounds
        psf['usepnt'] = self.usepnt
        psf['nxpix'] = nbins
        psf['nypix'] = nbins
        psf['binsz'] = wbin
        psf['coordsys'] = self.coord_sys
        psf['proj'] = self.proj
        if not self.usepnt:
            psf['xref'] = self.target[0] 
            psf['yref'] = self.target[1] 
        psf['amax'] = amax 
        psf['anumbins'] = abins
        if '.fits' in self.output:
            psf['logfile'] = self.output.replace('.fits', '.log')
        elif '.xml' in self.output:
            psf['logfile'] = self.output.replace('.xml', '.log')
        psf['debug'] = self.set_debug
        if self.set_log:
            psf.logFileOpen()
        if not self.__on_ram:
            psf.execute()
        else:
            psf.run()
        return

    # ctedispcube wrapper ---!
    def run_edispcube(self, cube, ebins=10, nbins=None, wbin=1.0, migramax=2.0, migrabins=100, ebin_alg='LOG', ebinfile=None, ebingamma=None, addbounds=False):
        '''Wrapper of ctedispcube.'''
        if nbins == None:
            nbins = int(self.roi * 2 / wbin)
        edisp = ctools.ctedispcube()
        edisp['inobs'] = self.input
        edisp['incube'] = cube
        edisp['caldb'] = self.caldb
        edisp['irf'] = self.irf
        edisp['outcube'] = self.output
        edisp['ebinalg'] = ebin_alg
        edisp['emin'] = self.e[0]
        edisp['emax'] = self.e[1]
        edisp['enumbins'] = ebins
        if ebinfile != None:
            edisp['ebinfile'] = ebinfile
        if ebingamma != None:
            edisp['ebingamma'] = ebingamma
        edisp['addbounds'] = addbounds
        edisp['usepnt'] = self.usepnt
        edisp['nxpix'] = nbins
        edisp['nypix'] = nbins
        edisp['binsz'] = wbin
        edisp['coordsys'] = self.coord_sys
        edisp['proj'] = self.proj
        if not self.usepnt:
            edisp['xref'] = self.target[0] 
            edisp['yref'] = self.target[1] 
        edisp['migramax'] = migramax 
        edisp['migrabins'] = migrabins
        if '.fits' in self.output:
            edisp['logfile'] = self.output.replace('.fits', '.log')
        elif '.xml' in self.output:
            edisp['logfile'] = self.output.replace('.xml', '.log')
        edisp['debug'] = self.set_debug
        if self.set_log:
            edisp.logFileOpen()
        if not self.__on_ram:
            edisp.execute()
        else:
            edisp.run()
        return

    # ctbkgcube wrapper ---!
    def run_bkgcube(self, cube, model):
        '''Wrapper of ctbkgcube.'''
        bkg = ctools.ctbkgcube()
        bkg['inobs'] = self.input
        bkg['inmodel'] = self.model
        bkg['incube'] = cube
        bkg['caldb'] = self.caldb
        bkg['irf'] = self.irf
        bkg['outcube'] = self.output
        bkg['outmodel'] = model    
        if '.fits' in self.output:
            bkg['logfile'] = self.output.replace('.fits', '.log')
        elif '.xml' in self.output:
            bkg['logfile'] = self.output.replace('.xml', '.log')
        bkg['debug'] = self.set_debug
        if self.set_log:
            bkg.logFileOpen()
        if not self.__on_ram:
            bkg.execute()       
        else: 
            bkg.run()
        return

    # ctlike wrapper ---!
    def run_maxlikelihood(self, binned=False, exp=None, bkg=None, psf=None, edisp=False, edispcube=None, fix_spat_for_ts=True):
        '''Wrapper of ctlike.'''
        if self.edisp:
            edisp = True
        like = ctools.ctlike()
        like['inobs'] = self.input
        like['inmodel'] = self.model
        if binned:
            like['expcube'] = exp
            like['psfcube'] = psf
            like['bkgcube'] = bkg
            if edisp:
                like['edispcube'] = edispcube
        if edisp:
            like['edisp'] = edisp
        like['outmodel'] = self.output
        like['caldb'] = self.caldb
        like['irf'] = self.irf
        like['refit'] = self.refit
        like['max_iter'] = self.max_iter
        like['like_accuracy'] = self.accuracy
        like['fix_spat_for_ts'] = fix_spat_for_ts
        like['statistic'] = self.stats
        like["nthreads"] = self.nthreads
        like['logfile'] = self.output.replace('.xml', '.log')
        like['debug'] = self.set_debug
        if self.set_log:
            like.logFileOpen()
        if not self.__on_ram: 
            like.execute()
        else:
            like.run()
        return

    # cterror wrapper ---!
    def run_asymerrors(self, asym_errors):
        '''Wrapper of cterror.'''
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
                err['debug'] = self.set_debug
                if self.set_log:
                    err.logFileOpen()
                if not self.__on_ram:
                    err.execute()
                else:
                    err.run()
        return self.output

    # ctulimit wrapper ---!
    def run_uplim(self):
        '''Wrapper of ctulimit.'''
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
        uplim['debug'] = self.set_debug
        if self.set_log:
            uplim.logFileOpen()
        if not self.__on_ram:
            uplim.execute()
        else:
            uplim.run()
        return

    # ctulimit wrapper ---!
    def run_lightcurve(self, nbins=20, bin_type='LIN'):
        '''Wrapper of cslightcrv.'''
        lc = cscripts.cslightcrv()
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
        lc['debug'] = self.set_debug
        if self.set_log:
            lc.logFileOpen()
        if not self.__on_ram:
            lc.execute()
        else:
            lc.run()
        return

    # cssens wrapper ---!
    def run_sensitivity(self, bins=20, wbin=0.05, enumbins=0):
        '''Wrapper of cssens.'''
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
        sens['debug'] = self.set_debug
        if self.set_log:
            sens.logFileOpen()
        if not self.__on_ram:
            sens.execute()
        else:
            sens.run()
        return

