# *****************************************************************************
# Copyright (C) 2021 INAF
# This software was provided as IKC to the Cherenkov Telescope Array Observatory
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
#
#    Simone Tampieri <simone.tampieri@studio.unibo.it>
#    Ambra Di Piano <ambra.dipiano@inaf.it>
#    Nicolò Parmiggiani <nicolo.parmiggiani@inaf.it>
#    Andrea Bulgarelli <andrea.bulgarelli@inaf.it>
#    Valentina Fioretti <valentina.fioretti@inaf.it>
#    Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
#    Antonio Addis <antonio.addis@inaf.it>
#    Giovanni De Cesare <giovanni.decesare@inaf.it>
#    Gabriele Panebianco <gabriele.panebianco3@unibo.it>
# *****************************************************************************

import numpy as np
from os import remove
from os.path import isfile
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from rtasci.aph import utils
from rtasci.aph.irf import aeff_eval
#from regions import CircleSkyRegion, Regions

# Bintable columns:
# 0  name = 'EVENT_ID'; format = '1J'; bscale = 1; bzero = 2147483648
#    name = 'TIME'; format = '1D'; unit = 's'
#    name = 'RA'; format = '1E'; unit = 'deg'
#    name = 'DEC'; format = '1E'; unit = 'deg'
#    name = 'ENERGY'; format = '1E'; unit = 'TeV'
# 5  name = 'DETX'; format = '1E'; unit = 'deg'
# 6  name = 'DETY'; format = '1E'; unit = 'deg'
#    name = 'MC_ID'; format = '1J'

class Photometrics():
    """This class contains the core of the RTAPH tool. It allows the extraction of on and off regions as well as perform aperture photometry."""
    def __init__(self, args):
        self.events_data = None
        self.events_filename = None
        self.mandatory_fields = ['RA', 'DEC', 'ENERGY']
        if 'events_filename' in args:
            self.events_filename = args['events_filename']
            self.events_data = self.load_data_from_fits_file(self.events_filename)
        elif 'events_list' in args:
            self.events_data = args['events_list']
        elif 'heatmap' in args:
            self.events_data = args['events_list']
        self.events_list_checks()

    def get_event_data_type(self):
        '''Return which event data is in use: observation list ("events_list") or events file ("events_filename").

        Parameter
        ---------

        Return
        ------
        data_input_type : str
            "events_list" for observation lists or "events_filename" for events files
        '''
        return type(self.events_data)

    def set_logger(self, logger):
        '''Setter for logger.
        
        Parameter
        ---------
        logger : logger obj
            logger for the class

        Return
        ------

        '''
        self.log = logger
        return self

    def events_list_checks(self):
        """Data con be a FITS_rec or a np.recarray
        see here: https://docs.astropy.org/en/stable/io/fits/usage/table.html

        Parameter
        ---------

        Return
        ------

        """
        if self.events_data is None:
            raise Exception('Events data is empy. Need a events list.')
        if isinstance(self.events_data, fits.fitsrec.FITS_rec):
            for f in self.mandatory_fields:
                if f not in self.events_data.columns.names:
                    raise Exception("Events data has no '{}' col".format(f))
        elif isinstance(self.events_data, np.recarray):
            for f in self.mandatory_fields:
                if f not in self.events_data.dtype.names:
                    raise Exception("Events data has no '{}' col".format(f))
        elif isinstance(self.events_data, np.ndarray):
            pass
        else:
            raise Exception("Events data must be FITS_rec or np.recarray")
        return self

    @staticmethod
    def load_data_from_fits_file(filename):
        """Load events extension data from a fits file.

        Parameter
        ---------
        filename : str
            path to a FITS file

        Returns
        -------
        data : ndarray 
            EVENTS data 
        """
        with fits.open(filename, mode='readonly') as hdul:
            data = hdul['EVENTS'].data
        return data

    def heatmap_region_counter(self, input_center, input_radius, binning, wcs=None):
        """Count photons in an input area.
        
        Parameter
        ---------
        input_center : dict or tuple or list
            region center coordinates with (ra, dec) in degrees
        input_radius : float or astropy Angle
            region radius in degrees
        
        Return
        ------
        counts : float
            counts in the region
        """
        region_center = utils.get_skycoord(input_center)
        region_radius = utils.get_angle(input_radius)

        # events coordinates from the selected events list
        ra, dec = np.arange(binning+1), np.arange(binning+1)
        if wcs is not None:
            ra, dec = wcs.world_to_pixel(SkyCoord(ra, dec, unit='deg', frame='icrs'))
        else:
            events_coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
        distances = region_center.separation(events_coords)
        return np.count_nonzero(distances < region_radius)

    def region_counter(self, input_center, input_radius, emin=None, emax=None, tmin=None, tmax=None):
        """Count photons in an input area.
        
        Parameter
        ---------
        input_center : dict or tuple or list
            region center coordinates with (ra, dec) in degrees
        input_radius : float or astropy Angle
            region radius in degrees
        emin : float or None
            minimum energy selection cut
        emax : float or None
            maximum energy selection cut
        tmin : float or None
            minimum time selection cut
        tmax : float or None
            maximum time selection cut
        
        Return
        ------
        counts : float
            counts in the region
        """
        region_center = utils.get_skycoord(input_center)
        region_radius = utils.get_angle(input_radius)

        # filtering...
        condlist = np.full(len(self.events_data.field('ENERGY')), True)
        # ... w/ energy boundaries
        if emin is not None:
            condlist &= self.events_data.field('ENERGY') >= emin
        if emax is not None:
            condlist &= self.events_data.field('ENERGY') <= emax
        # FIXME: TIME needs a better implementation
        # atm it consider users that knows the time format in the input fits
        if tmin is not None:
            condlist &= self.events_data.field('TIME') >= tmin
        if tmax is not None:
            condlist &= self.events_data.field('TIME') <= tmax

        events_list = np.extract(condlist, self.events_data)
        # events coordinates from the selected events list
        events_coords = SkyCoord(events_list.field('RA'), events_list.field('DEC'), unit='deg', frame='icrs')
        distances = region_center.separation(events_coords)
        return np.count_nonzero(distances < region_radius)

    @classmethod
    def coords_to_absolute(self, coords):
        '''Take coordinates to absolute value.

        Parameter
        ---------
        coords : dict
            coordinates with ra, dec in degrees

        Return
        ------
        coords : dict
            coordinates with ra, dec absolute value
        '''
        if np.sign(coords['ra']) < 0:
            self.ra_neg = True
            coords['ra'] = np.abs(coords['ra'])
        else:
            self.ra_neg = False
        if np.sign(coords['dec']) < 0:
            self.dec_neg = True
            coords['dec'] = np.abs(coords['dec'])
        else:
            self.dec_neg = False
        return coords

    @classmethod
    def regions_reassign_sign(self, regions):
        '''Reassign sign to coordinates.
        
        Parameter
        ---------
        regions : list
            list of regions dictionaries with ra, dec in ansolute value

        Return
        ------
        regions : list
            list of regions dictionaries with ra, dec in degrees
        '''
        for r in regions:
            if self.ra_neg:
                r['ra'] = -r['ra']
            if self.dec_neg:
                r['dec'] = -r['dec']
        return regions

    @classmethod
    def reflected_regions_with_astropy(self, input_pointing_center, input_region_center, input_region_radius):
        """Find regions with reflected algorithm with astropy projection.

        Parameter
        ---------
        input_pointing_center : astropy SkyCoord or dict
            pointing coordinates in degrees
        input_region_center : astropy SkyCoord or dict
            on region center coordinates in degrees
        input_region_radius : astropy Angle or float
            region radius in degrees

        Returns
        -------
        regions : list
            list of regions dictionary with ra, dec, rad in degrees
        """
        pointing_center = utils.get_skycoord(input_pointing_center)
        region_center = utils.get_skycoord(input_region_center)
        region_radius = utils.get_angle(input_region_radius)

        # Angular separation of reflected regions. 1.05 factor is to have a margin
        region_diameter = 1.05 * 2.0 * region_radius
        radius = pointing_center.separation(region_center)
        # the numbers_of_reflected regions is the number of center that can stay
        # on the circumference, NOT the really computated number of regions.
        numbers_of_reflected_regions = int(2 * np.pi * radius / region_diameter)
        # Indeed, we skip the source region and the two near (see below), so we
        # need at least 4 centers to get one off region.
        if numbers_of_reflected_regions < 4:
            raise Exception('the combination of region radius and coordinates does not allow to compute reflected regions.')
        regions_offset_angle = Angle(360, unit='deg') / numbers_of_reflected_regions
        
        regions = []
        # starting from the source region 0, we skip region 1 and region N, so 2..N-1
        starting_pos_angle = pointing_center.position_angle(region_center)
        for i in range(2, numbers_of_reflected_regions-1):
            theta = starting_pos_angle + i * regions_offset_angle
            coord_pos = pointing_center.directional_offset_by(theta, radius)
            regions.append({ 'ra': coord_pos.ra.deg, 'dec': coord_pos.dec.deg, 'rad': region_radius.deg })
        return regions

    @classmethod
    def reflected_regions(self, pointing, target, region_radius, skip_adjacent=True, min_regions_number=4):
        """Find regions with reflected algorithm.
        
        Parameter
        ---------
        pointing : astropy SkyCoord or dict
            pointing coordinates in degrees
        target : astropy SkyCoord or dict
            target coordinates in degrees
        region radius : astropy Angle or float
            region radius in degrees

        Returns
        -------
        regions : list
            list of regions dictionary with ra, dec, rad in degrees
        """

        pointing = utils.get_skycoord(pointing)
        target = utils.get_skycoord(target)        
        offaxis_angle = utils.get_offset(pointing=(pointing.ra.deg, pointing.dec.deg), target=(target.ra.deg, target.dec.deg))

        # Angular separation of reflected regions. 1.05 factor is to have a margin
        region_diameter = 1.05 * 2.0 * region_radius
        numbers_of_reflected_regions = int(2 * np.pi * offaxis_angle / region_diameter)
        
        # Skip the source region and the two near => at least 4 centers to get one off region.
        if numbers_of_reflected_regions < min_regions_number:
            raise Exception('the combination of region radius and coordinates does not allow to compute reflected regions.')
        regions_offset_angle = 360 / numbers_of_reflected_regions
        
        regions = []
        # starting from the source region 0, we skip region 1 and region N, so 2..N-1
        starting_pos_angle = utils.angle_between(pointing=(pointing.ra.deg, pointing.dec.deg), target=(target.ra.deg, target.dec.deg))
        if skip_adjacent:
            start_reg_count = 2
            stop_reg_count = numbers_of_reflected_regions-1
        else:
            start_reg_count = 1
            stop_reg_count = numbers_of_reflected_regions
        for i in range(start_reg_count, stop_reg_count):
            theta = starting_pos_angle + i * regions_offset_angle
            ra = offaxis_angle * np.cos(np.deg2rad(theta)) + pointing.ra.deg
            dec = offaxis_angle * np.sin(np.deg2rad(theta)) + pointing.dec.deg
            regions.append({ 'ra': ra, 'dec': dec, 'rad': region_radius })
        return regions

    @classmethod
    def wobble_regions_with_astropy(self, *args):
        '''Find wobble regions as cross regions with astropy projection.

        Parameter
        ---------

        Return
        ------
        regions : list
            list of region dictionaries with ra, dec, rad in degrees
        '''
        return self.cross_regions_with_astropy(*args)

    @classmethod
    def wobble_regions(self, *args):
        '''Find wobble regions as cross regions.

        Parameter
        ---------

        Return
        ------
        regions : list
            list of region dictionaries with ra, dec, rad in degrees
        '''
        return self.cross_regions(*args)

    @classmethod
    def cross_regions_with_astropy(self, input_pointing_center, input_region_center, input_region_radius):
        """Return the three background regions starting from pointing and source one with astropy projection.

        Parameters
        ----------
        input_pointing_center : astropy SkyCoord or dict
            pointing coordinates in degrees
        input_region_center : astropy SkyCoord or dict
            on region coordinates in degrees
        input_region_radius : astropy Angle or float
            region radius in degrees

        Returns
        -------
        regions : list
            list of regions dictionary with ra, dec, rad in degrees
        """
        # FIXME Wobble algorithm has no check about distance and region radius.
        pointing_center = utils.get_skycoord(input_pointing_center)
        region_center = utils.get_skycoord(input_region_center)
        region_radius = utils.get_angle(input_region_radius)
        radius = pointing_center.separation(region_center)
        starting_pos_angle = pointing_center.position_angle(region_center)
        regions = []
        for i in range(1,4):
            theta = starting_pos_angle + i * Angle(90, unit='deg')
            coord_pos = pointing_center.directional_offset_by(theta, radius)
            regions.append({ 'ra': coord_pos.ra.deg, 'dec': coord_pos.dec.deg, 'rad': region_radius.deg })
        return regions

    @classmethod
    def cross_regions(self, pointing, target, region_radius):
        """Return the three background regions starting from pointing and source one.

        Parameters
        ----------
        pointing : astropy SkyCoord or dict
            pointing coordinates in degrees
        target : astropy SkyCoord or dict
            target coordinates in degrees
        region_radius : astropy Angle or float
            region radius in degrees

        Returns
        -------
        regions : list
            list of regions dictionary with ra, dec, rad in degrees
        """
        pointing = utils.get_skycoord(pointing)
        target = utils.get_skycoord(target)        
        offaxis_angle = utils.get_offset(pointing=(pointing.ra.deg, pointing.dec.deg), target=(target.ra.deg, target.dec.deg))
        
        # starting from the source region 0, we skip region 1 and region N, so 2..N-1
        starting_pos_angle = utils.angle_between(pointing=(pointing.ra.deg, pointing.dec.deg), target=(target.ra.deg, target.dec.deg))
        regions = []
        for i in range(1,4):
            theta = starting_pos_angle + i * 90
            ra = offaxis_angle * np.cos(np.deg2rad(theta)) + pointing.ra.deg
            dec = offaxis_angle * np.sin(np.deg2rad(theta)) + pointing.dec.deg
            regions.append({ 'ra': ra, 'dec': dec, 'rad': region_radius })
        return regions

    @classmethod
    def write_region(self, coords, filename, **kwargs):
        '''Write region file.

        Parameter
        ---------
        coords : list
            list of regions dictionary coordinates
        filename : str
            path to destination file

        Return
        ------

        '''
        try:
            iter(coords)
        except TypeError:
            raise Exception('Coords must be iterable')
        circles = []
        for coord in coords:
            center = utils.get_skycoord(coord)
            rad = utils.get_angle(coord['rad'])
            circles.append(CircleSkyRegion(center=center, radius=rad, visual=kwargs))
        off_regions = Regions(circles) 
        off_regions.serialize(format='ds9')
        if isfile(filename): 
            remove(filename)
        off_regions.write(filename, format='ds9')
        return self


    def counting(self, src, rad, off_regions, e_min=None, e_max=None, t_min=None, t_max=None, draconian=False):
        '''Count all photometric results.
        
        Parameter
        ---------
        src : dict
            source coordinates ra,dec in degrees
        rad : float
            region radius in degrees
        off_regions : list
            list of regions dictionary with ra, dec, rad in degrees
        e_min : float or None
            minimum energy selection cut
        e_max : float or None
            maximum energy selection cut
        t_min : float or None
            minimum time selection cut
        t_max : float or None
            maximum time selection cut
        draconian : bool
            instead of returning error note raise error if on/off counts < 10

        Return
        ------
        on_counts : float
            source counts
        off_counts : float
            background counts
        alpha : float
            alpha parameter
        excess : float
            excess counts
        signif : float
            Li & Ma significance in gaussian sigmas
        err_note : str or None
            error message
        '''
        on_count = self.region_counter(src, rad, emin=e_min, emax=e_max, tmin=t_min, tmax=t_max)
        off_count = 0
        # compute off counts
        for r in off_regions:
            off_count += self.region_counter(r, r['rad'], emin=e_min, emax=e_max, tmin=t_min, tmax=t_max)

        alpha = 1 / len(off_regions)
        excess = on_count - alpha * off_count
        signif = utils.li_ma(on_count, off_count, alpha)
        # !!! here we can implement checks
        err_note = None
        # handle error 
        if on_count < 10 or off_count < 10:
            err_note = 'Not compliant with Li & Ma requirements.'
            if draconian:
                raise Exception(err_note)
        return on_count, off_count, alpha, excess, signif, err_note

    def find_off_regions(self, algo, src, pnt, rad, skip_adjacent=True, verbose=False, save=False):
        '''Produces off regions for aperture photometry.
        
        Parameter
        ---------
        algo : str
            algorithm for off regions determination
        src : dict
            source coordinates with ra, dec in degrees
        pnt : dict
            pointing coordinates with ra, dec in degrees
        rad : float
            region radius in degrees
        verbose : bool
            message verbosity
        save : bool
            save regions to file
        '''
        if not len(pnt) == 2 and len(src) == 2:
            raise Exception('need source and pointing coordinates and a region radius to do aperture photometry')

        off_regions = None
        if algo.lower() == 'cross' or algo.lower() == 'wobble':
            off_regions = self.cross_regions(pnt, src, rad)
        elif algo.lower() == 'reflection':
            off_regions = self.reflected_regions_with_astropy(pnt, src, rad)
        else:
            raise Exception('invalid background regions algorithm')

        if save:
            self.write_region(off_regions, save, color='red', dash=True, width=2)
        return off_regions

class PhmConfiguration(object):
    """This class allows to obtain a compatible configuration object for Photometrics class.
    
    Parameter
    ---------
    d : dict
        configuration paramers key-value pairs
    """
    def __init__(self, d):
        self.__dict__ = d
        self.__dict__['begin_time'] = float(self.__dict__['begin_time'])
        self.__dict__['end_time'] = float(self.__dict__['end_time'])
        self.__dict__['source_ra'] = float(self.__dict__['source_ra'])
        self.__dict__['source_dec'] = float(self.__dict__['source_dec'])
        self.__dict__['region_radius'] = float(self.__dict__['region_radius'])
        self.__dict__['verbose'] = float(self.__dict__['verbose'])
        self.__dict__['emin'] = float(self.__dict__['emin'])
        self.__dict__['emax'] = float(self.__dict__['emax'])
        self.__dict__['pixel_size'] = float(self.__dict__['pixel_size'])
        self.__dict__['power_law_index'] = float(self.__dict__['power_law_index'])
        self.__dict__['irf_file'] = str(self.__dict__['irf_file'])
        
def photometrics_counts(events_list, events_type, pointing, target_coords, region_rad=0.2, skip_adjacent=True):
    '''Return photometric counts.
    
    Parameter
    ---------
    events_list : str
        path to events list FITS file
    events_type : str   
        observation list or events file
    pointing : dict
        pointing coordinates with ra, dec in degrees
    target_coords : dict
        on region coordinates with ra, dec in degrees

    Return
    ------
    photometry : dict
        photometric results with on, off, excess counts and alpha parameter
    '''
    phm = Photometrics({events_type: events_list})
    reflected_regions = phm.reflected_regions(pointing, target_coords, region_rad, skip_adjacent)
    on_count = phm.region_counter(target_coords, region_rad)
    off_count = 0
    for r in reflected_regions:
        off_count += phm.region_counter(r, r['rad'])
    alpha = 1 / len(reflected_regions)
    return {'on': on_count, 'off': off_count, 'alpha': alpha, 'excess': on_count - alpha * off_count}

def phm_options(erange, trange, target, pointing, irf_file, index=-2.4, bkg_method="reflection", pixsize=0.05, verbose=0, save_off_reg='.'):
    '''Return arguments dictionary for Photometrics class.
    
    Parameter
    ---------
    erange : list
        energy interval with [emin, emax] in TeV
    trange : list
        time interval with [tmin, tmax] in s
    target : dict
        target dictionary with ra, dec, rad in degrees
    pointing : dict
        pointing dictionary with ra, dec in degrees
    irf_file : str
        path to IRF file in FITS format
    index : float
        spectral index of powerlaw
    bkg_methos : str
        aperture photometry bkg estimation method
    pixsize : float
        size of a pixe in degrees
    verbose : int
        verbose quantity
    save_off_region : str
        path to destination file

    Return
    ------
    configuration : obj Phm
        Photometrics configuration
    '''
    opts = {}
    opts['verbose'] = verbose
    opts['irf_file'] = irf_file
    opts['source_ra'] = target['ra']
    opts['source_dec'] = target['dec']
    opts['pointing_ra'] = pointing['ra']
    opts['pointing_dec'] = pointing['dec']
    opts['region_radius'] = target['rad']
    opts['background_method'] = bkg_method
    opts['save_off_regions'] = save_off_reg
    opts['emin'] = erange[0]
    opts['emax'] = erange[1]
    opts['pixel_size'] = pixsize
    opts['begin_time'] = trange[0]
    opts['end_time'] = trange[1]
    opts['step_time'] = trange[1]-trange[0]
    opts['power_law_index'] = index
    configuration = PhmConfiguration(opts)
    return configuration


def get_counts_in_region(filename, target, pointing, trange, erange, off_regions='off_regions.reg'):
    '''Get on counts in target region.
    
    Parameter
    ---------
    filename : str
        path to events list FITS file
    target : dict
        target coordinates with ra, dec, rad in degrees
    pointing : dict
        pointing coordinates with ra, dec in degrees
    trange : list
        time interval with [tmin,temax] in s
    erange : list
        energy interval with [emin, emax] in TeV
    off_regions : str
        path to destination file
    
    Return
    ------
    on :  float
        on counts
    '''
    phm = Photometrics({'events_filename': filename})
    off_regions = phm.find_off_regions(algo='cross', src=target, pnt=pointing, rad=target['rad'], save=off_regions)
    on, off, alpha, excess, sigma, err_note = phm.counting(src=target, rad=target['rad'], off_regions=off_regions, e_min=erange[0], e_max=erange[1], t_min=trange[0], t_max=trange[1], draconian=False)
    return on

def get_excess_in_region(filename, target, pointing, trange, erange, off_regions='off_regions.reg'):
    '''Get excess counts of target region.
    
    Parameter
    ---------
    filename : str
        path to events list FITS file
    target : dict
        target coordinates with ra, dec, rad in degrees
    pointing : dict
        target coordinates with ra, dec in degrees
    trange : list
        time interval [tmin, tmax] in s
    erange : list
        energy interval [emin, emax] in TeV
    off_regions : str
        path to destination file

    Return
    ------
    excess : float
        excess counts
    '''
    phm = Photometrics({'events_filename': filename})
    off_regions = phm.find_off_regions(algo='cross', src=target, pnt=pointing, rad=target['rad'], save=off_regions)
    on, off, alpha, excess, sigma, err_note = phm.counting(src=target, rad=target['rad'], off_regions=off_regions, e_min=erange[0], e_max=erange[1], t_min=trange[0], t_max=trange[1], draconian=False)
    return excess

def get_off_counts_wrt_region(filename, target, pointing, trange, erange, off_regions='off_regions.reg'):
    '''Get off counts in target region.
    
    Parameter
    ---------
    filename : str
        path to events list FITS file
    target : dict
        target coordinates with ra, dec, rad in degrees
    pointing : dict
        target coordinates with ra, dec in degrees
    trange : list
        time interval [tmin, tmax] in s
    erange : list
        energy interval [emin, emax] in TeV
    off_regions : str
        path to destination file

    Return
    ------
    off : float
        off counts
    '''
    phm = Photometrics({'events_filename': filename})
    off_regions = phm.find_off_regions(algo='cross', src=target, pnt=pointing, rad=target['rad'], save=off_regions)
    on, off, alpha, excess, sigma, err_note = phm.counting(src=target, rad=target['rad'], off_regions=off_regions, e_min=erange[0], e_max=erange[1], t_min=trange[0], t_max=trange[1], draconian=False)
    return off

def get_exposure_in_region(target, pointing, trange, erange, irf, index, bkg_method='cross'):
    '''Get exposure (cm2/s) in target region.
    
    Parameter
    ---------
    target : dict
        target coordinates with ra, dec, rad in degrees
    pointing : dict
        pointing coordinates with ra, dec in degrees
    trange : list
        time interval with [tmin, tmax] in s
    erange : list
        energy interval with [emin, emax] in TeV
    irf : str
        path to IRF file in FITS format
    index : float
        spectral index of powerlaw
    bkg_method : str
        aperture photometry background estimation method

    Return
    ------
    exposure : float
        region exposure in cm2/s
    '''
    args = phm_options(erange=erange, trange=trange, target=target, pointing=pointing, irf_file=irf, index=index, bkg_method=bkg_method)
    region_eff_resp = aeff_eval(args, src=target, pnt=pointing)
    livetime = trange[1]-trange[0]
    exposure = region_eff_resp * livetime
    return exposure

def get_aeff_in_region(target, pointing, trange, erange, irf, index, bkg_method='cross'):
    '''Get effective area (cm2) in target region.
 
    Parameter
    ---------
    target : dict
        target coordinates with ra, dec, rad in degrees
    pointing : dict
        pointing coordinates with ra, dec in degrees
    trange : list
        time interval with [tmin, tmax] in s
    erange : list
        energy interval with [emin, emax] in TeV
    irf : str
        path to IRF file in FITS format
    index : float
        spectral index of powerlaw
    bkg_method : str
        aperture photometry background estimation method

    Return
    ------
    exposure : float
        effective area in cm2
    '''
    args = phm_options(erange=erange, trange=trange, target=target, pointing=pointing, irf_file=irf, index=index, bkg_method=bkg_method)
    region_eff_resp = aeff_eval(args, src=target, pnt=pointing)
    return region_eff_resp



def get_prefactor_from_bkg_and_sigma(sigma, bkg, target, pointing, trange, erange, irf, spectral_index):
    '''From sigma and bkg find the required prefactor.
    
    Parameter
    ---------
    sigma : float
        signal significance
    bkg : float
        background counts
    target : dict
        target coordinates with ra, dec, rad in degrees
    pointing : dict
        pointing coordinates with ra, dec in degrees
    trange : list
        time interval with [tmin, tmax] in s
    erange : list
        energy interval with [emin, emax] in TeV
    irf : str
        path to IRF file in FIST format
    spectral_index : float
        spectral index of powerlaw

    Return
    ------
    prefactor : float
        prefactor / amplitude / normalisation of powerlaw in ph/cm2/s/MeV
    '''
    excess = utils.get_excess_from_sigma_and_bkg(sigma=sigma, bkg=bkg)
    exposure = get_exposure_in_region(target=target, pointing=pointing, trange=trange, erange=erange, irf=irf, index=spectral_index)
    flux = excess / exposure 
    prefactor = utils.get_prefactor(flux=flux, erange=erange, gamma=spectral_index, unit='TeV')
    return prefactor
