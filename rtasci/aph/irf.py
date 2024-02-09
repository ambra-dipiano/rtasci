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

import math
import numpy as np
from scipy import interpolate, integrate
from rtasci.aph import utils 
from astropy.coordinates import SkyCoord
from astropy.io import fits

class IRF:
    """Class that allows to operate on the Instrument Response Functions in FITS format.
    
    Parameter
    ---------
    filename : str
        path to IRF in FITS format
    """
    def __init__(self, filename):
        self.filename = filename
        self.hdul = fits.open(self.filename)

    def get_extension(self, name):
        """Extract extension content by keyword.
        
        Parameter
        ---------
        name : str
            keyword of the FITS extension
        """
        return self.hdul[name]

    def get_eff_area(self):
        """Extract the effective area extension.

        Parameter
        ---------

        Return
        ------
        aeff : FITS extension
            effective area FITS extension content
        """
        return self.get_extension('EFFECTIVE AREA')

    def get_psf_data(self):
        """Extract the point spread function extension.
        
        Parameter
        ---------

        Return
        ------
        psf : FITS extension
            point spread function FITS extension content
        """
        return self.get_extension('POINT SPREAD FUNCTION')

class EffectiveArea:
    """Class that specifically operates on the Effective Area extension of an Instrument Response Function in FITS format.
    
    Parameter
    ---------
    irf_filename : str
        path to the IRF file in FITS format
    eff_area_bintable : FITS extension
        effective area FITS extension
    """
    def __init__(self, irf_filename=None, eff_area_bintable=None):
        self.irf_filename = None
        self.eff_area = None
        # a sort of cache...
        self.energies = None
        self.thetas = None
        self.aeff_matrix = None
        # check path
        if irf_filename is not None:
            self.irf_filename = irf_filename
            irf = IRF(self.irf_filename)
            self.eff_area = irf.get_eff_area()
        # check extension
        elif eff_area_bintable is not None:
            self.eff_area = eff_area_bintable
        # check table
        if self.eff_area is None:
            raise Exception('Need an irf or effective area bintable')
        # get data
        self.get_data_matrices()

    def columns(self):
        """Get FITS extension columns.
        
        Parameter
        ---------
        
        Return
        ------
        columns : list
            FITS extension colums
        """
        return self.eff_area.columns

    def data(self):
        """Get FITS extension data.
        
        Parameter
        ---------

        Return
        ------
        data : ndarray
            FITS extension data
        """
        return self.eff_area.data

    def get_data_matrices(self):
        """Get effective area data matrix, energy[LO,HI] and theta[LO,HI].

        Parameter
        ---------

        Return
        ------
        aeff_matrix : ndarray
            effective area FITS extension matrix
        energies : ndarray
            energy FITS extension subset
        thetas : ndarray
            thetas FITS extension subset
        """
        if self.energies is None or self.thetas is None or self.aeff_matrix is None:
            data = self.data()
        if self.aeff_matrix is None:
            self.aeff_matrix = data.field('EFFAREA')[0]
        if self.energies is None:
            self.energies = np.column_stack((data.field('ENERG_LO')[0], data.field('ENERG_HI')[0]))
        if self.thetas is None:
            self.thetas = np.column_stack((data.field('THETA_LO')[0], data.field('THETA_HI')[0]))
        return self.aeff_matrix, self.energies, self.thetas

    def get_aeff_1d_log(self, offset, energy):
        """Compute the effective area in [m2] via 1d interpolation on energy.

        Parameters
        ----------
        offset : astropy Angle
            region offset from the pointing in astropy Angle
        energy : float
            reference energy in [TeV]

        Return
        ------
        aeff : float
            effective area value in [m2]
        """
        offset_angle = utils.get_angle(offset)
        aeff_matrix, energy_bins, theta_bins = self.get_data_matrices()

        theta_index = None
        for i, tb in enumerate(theta_bins):
            if offset_angle.degree >= tb[0] and offset_angle.degree < tb[1]:
                theta_index = i
                break
        if theta_index is None:
            raise Exception('Theta offset is out of range ({})'.format(offset))

        # energy interpolation
        energy_mid = [ (e[0]+e[1])/2 for e in np.log10(energy_bins) ]
        energy_fn  = interpolate.interp1d(x = energy_mid, y = aeff_matrix[theta_index])
        return energy_fn(np.log10(energy))

    def get_aeff_2d_log(self, input_offset, input_energy):
        """Compute the effective area in [m²] via 2d interpolation on thetas and energy.

        Parameters
        ----------
        offset : astropy Angle
            region offset from the pointing in astropy Angle
        energy : float
            reference energy in [TeV]

        Return
        ------
        aeff : float
            effective area value in [m2]
        """
        offset_angle = utils.get_angle(input_offset)
        aeff_matrix, energy_bins, theta_bins = self.get_data_matrices()

        # energy interpolation
        theta_mid  = [ (t[0]+t[1])/2 for t in theta_bins ]
        energy_mid = [ (e[0]+e[1])/2 for e in np.log10(energy_bins) ]
        energy_fn  = interpolate.interp2d(x = energy_mid, y = theta_mid, z = aeff_matrix)
        return energy_fn(np.log10(input_energy), input_offset)[0]

    def weighted_value_for_region(self, *args):
        """Compute the effective response value weighter over the photometric region.
        
        Parameter
        ---------
        *args

        Return
        ------
        aeff : float
            effective area value in [m2]
        """
        return self.weighted_aeff_flat_psf_w_powerlaw(*args)

    def weighted_aeff_flat_psf_w_powerlaw(self, region, pointing, input_energies, pixel_size=0.05, e_index=-2.4):
        """Compute effective area value [m²] for a specific region. The idea is integrate the PSF gradually but the offset from center and from pointing works differently we need to associate the PSF degradation from source region center and the AEFF degradation form the pointing. This method uses the PSF value plain as it is easy for energy weight but it does not consider the spatial component of the PSF.

        Parameters
        ----------
        region : dict
            region dictionary with ra, dec, rad in degrees
        pointing : dict
            pointing dictionary with ra, dec in degrees
        energies : list
            energy interval in [TeV] given as [emin, emax]
        pixel_size : float
            a pixel size in degrees
        e_index : float
            powerlaw index value

        Return
        ------
        val : float
            effective area value in [m2]
        """
        if len(input_energies) != 2:
            raise Exception('need two energies')

        psf = PSF(irf_filename=self.irf_filename)

        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        log_energies = np.log10(input_energies)
        # N steps for every unit of log energy
        steps = int(np.ceil(log_energies[1]-log_energies[0]) * 10)
        energies = 10**np.linspace(log_energies[0], log_energies[1], steps)
        powerlaw = lambda x: x**e_index
        i_full = integrate.quad(powerlaw, input_energies[0], input_energies[1])
        i_partials = [ integrate.quad(powerlaw, energies[i], energies[i+1]) for i,v in enumerate(energies[:-1]) ]
        i_factor = [ p[0]/i_full[0] for p in i_partials ]
        energies_middle = (energies[1:]+energies[:-1])/2

        psf_engines = {}
        for en in energies_middle:
            psf_engines[en] = psf.get_psf_engine(region, pointing, en)
        region_radius_rad = np.deg2rad(region['rad'])

        n_points = len(offsets)
        val = 0
        for t in offsets:
            for i, en in enumerate(energies_middle):
                psf_rate = psf_engines[en](0, region_radius_rad)[0]
                val += self.get_aeff_2d_log(t, en) * i_factor[i] * psf_rate / n_points
        return val

    def weighted_value_for_region_w_powerlaw(self, region, pointing, input_energies, pixel_size=0.05, e_index=-2.4):
        """Compute the effective area value in [m²] for a specific region. This method uses the energy range to evaluate the AEFF. The energy range is binned and weighted with a powerlaw with index. Each pixel has a radial weigth and a specific column of energies weight. Lower energies have more weigth than higher. If the energy range is small, the effect is trascurable similarly to weighted_value_for_region_single_energy method.

        Parameters
        ----------
        region : dict
            region dictionary with ra, dec, rad in degrees
        pointing : dict
            pointing dictionary with ra, dec in degrees
        energies : list
            energy interval in [TeV] given as [emin, emax]
        pixel_size : float
            size of a pixel in degrees
        e_index : float
            powerlaw spectral index

        Return
        ------
        val : float
            effective area value in [m2] 
        """
        if len(input_energies) != 2:
            raise Exception('need two energies')

        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        log_energies = np.log10(input_energies)
        # N steps for every unit of log energy
        steps = int(np.ceil(log_energies[1]-log_energies[0]) * 10)
        energies = 10**np.linspace(log_energies[0], log_energies[1], steps)
        powerlaw = lambda x: x**e_index
        i_full = integrate.quad(powerlaw, input_energies[0], input_energies[1])
        i_partials = [ integrate.quad(powerlaw, energies[i], energies[i+1]) for i,v in enumerate(energies[:-1]) ]
        i_factor = [ p[0]/i_full[0] for p in i_partials ]
        energies_middle = (energies[1:]+energies[:-1])/2
        n_points = len(offsets)
        val = 0
        for t in offsets:
            for i, en in enumerate(energies_middle):
                val += self.get_aeff_2d_log(t, en) * i_factor[i] / n_points
        return val

    def weighted_value_for_region_no_powerlaw(self, region, pointing, input_energies, pixel_size=0.05):
        """Compute the effective area value [m²] for a specific region. This method use an energy range to evaluate the aeff. The energy range is binned and every matrix cube (pixel distance * energy bin) have the same weight. No powerlaw is considered. 

        Parameters
        ----------
        region : dict
            region dictionary with ra, dec, rad in degrees
        pointing : dict
            pointing dictionary with ra, dec in degrees
        energies : list
            energy interval in [TeV] given as [emin, emax]
        pixel_size : float
            size of a pixel in degrees

        Return
        ------
        val : float
            effective area value in [m2]
        """
        if len(input_energies) != 2:
            raise Exception('need two energies')

        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        log_energies = np.log10(input_energies)
        diff = np.ceil(log_energies[1]-log_energies[0])
        steps = int(diff * 10) # N steps for every unit of log energy
        energies = 10**np.linspace(log_energies[0], log_energies[1], steps)
        n_points = len(offsets) * len(energies)
        val = 0
        for t in offsets:
            for en in energies:
                val += self.get_aeff_2d_log(t, en) / n_points
        return val

    def weighted_value_for_region_single_energy(self, region, pointing, energy, pixel_size=0.05):
        """Compute the effective area value [m²] for a specific region. This method is just a plain output from an array of points and one energy (usually the middle point of the energy range. DEPRECATED.

        Parameters
        ----------
        region : dict
            region dictionary with ra, dec, rad in degrees
        pointing : dict
            pointing dictionary with ra, dec in degrees
        energies : list
            energy interval in [TeV] given as [emin, emax]
        pixel_size : float
            size of a pixel in degrees

        Return
        ------
        val : float
            effective area value in [m2]
        """
        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        n_points = len(offsets)
        val = 0
        for t in offsets:
            val += self.get_aeff_2d_log(t, energy) / n_points
        return val # m2

    # helpers
    @staticmethod
    def create_pixel_map(region, pixel_side):
        """Compute the map pixel grid shifting from side to center.

        Paramter
        --------
        pixel_side : array
            edges of pixels

        Return
        ------
        pixels_midpoint : list
            centers of pixels
        """
        for k in ['ra', 'dec', 'rad']:
            if k in region:
                continue
            raise Exception('region data missing {} mandatory key.'.format(k))
        if region['rad'] <= 0:
            raise Exception('region radius must be > 0')
        if pixel_side <= 0:
            raise Exception('pixel side must be > 0')

        region_center = { 'ra': float(region['ra']), 'dec': float(region['dec']) }
        region_rad = utils.get_angle(float(region['rad']))
        pixel_side_angle = utils.get_angle(float(pixel_side))

        # +10% to get a bit of margin
        n_pixel_on_rad = 1.1* region_rad / pixel_side_angle
        if n_pixel_on_rad <= 1:
            n_pixel_on_rad = 1

        n_pixel_on_axis = float(math.ceil(n_pixel_on_rad))
        multipliers = np.arange(-1*n_pixel_on_axis, n_pixel_on_axis+1)
        pixels_midpoint = []
        for i in multipliers:
            for j in multipliers:
                pixels_midpoint.append({ 'ra':  region_center['ra']  + i * pixel_side_angle.deg,
                                         'dec': region_center['dec'] + j * pixel_side_angle.deg })
        return pixels_midpoint

    @staticmethod
    def select_points_in_region(midpoints, region):
        """Select pixels within region.

        Parameter
        ---------
        midpoints : array
            center point of pixels
        region : dict
            region dictionary with ra, dec, rad in degrees

        Return
        ------
        midpoints : array
            selected array of midpoints within region
        """
        for k in ['ra', 'dec', 'rad']:
            if k in region:
                continue
            raise Exception('region data missing {} mandatory key.'.format(k))
        if region['rad'] <= 0:
            raise Exception('region radius must be > 0')
        if len(midpoints) < 1:
            raise Exception('need at least 1 point to check')

        region_center = utils.get_skycoord(region)
        region_radius = utils.get_angle(region['rad'])
        midpoints_ra = []
        midpoints_dec = []
        for p in midpoints:
            midpoints_ra.append(p['ra'])
            midpoints_dec.append(p['dec'])
        midpoints_coords = SkyCoord(midpoints_ra, midpoints_dec, unit='deg', frame='icrs')
        distances = region_center.separation(midpoints_coords)
        return np.extract(distances < region_radius, midpoints)

    @staticmethod
    def get_thetas(point, midpoints):
        """Get thetas of points in region.

        Parameter
        ---------
        point : dict
            center of region with ra, dec in degrees
        midpoints : array
            pixel centers

        Return
        ------
        thetas : list
            separation of midpoints from region center
        """
        for k in ['ra', 'dec']:
            if k in point:
                continue
            raise Exception('point coord {} is missing.'.format(k))
        if len(midpoints) < 1:
            raise Exception('need at least 1 point to check')
        pnt = utils.get_skycoord(point)
        midpoints_ra = []
        midpoints_dec = []
        for p in midpoints:
            midpoints_ra.append(p['ra'])
            midpoints_dec.append(p['dec'])
        midpoints_coords = SkyCoord(midpoints_ra, midpoints_dec, unit='deg', frame='icrs')
        return [ ang.degree for ang in pnt.separation(midpoints_coords) ]
    
class PSF:
    """Class that specifically operates on the Effective Area extension of an Instrument Response Function in FITS format.
    
    Parameter
    ---------
    irf_filename : str
        path to the IRF file in FITS format
    psf_area_bintable : FITS extension
        point spread function FITS extension
    """
    def __init__(self, irf_filename=None, psf_bintable=None):
        self.irf_filename = None
        self.psf_data = None
        self.fields = ('SIGMA_1', 'SIGMA_2', 'SIGMA_3', 'SCALE', 'AMPL_2', 'AMPL_3')
        # a sort of cache...
        self.energies = None
        self.thetas = None
        self.psf_matrix = None
        # check file
        if irf_filename is not None:
            self.irf_filename = irf_filename
            irf = IRF(self.irf_filename)
            self.psf_data = irf.get_psf_data()
        # check table
        elif psf_bintable is not None:
            self.psf_data = psf_bintable
        # check extension data
        if self.psf_data is None:
            raise Exception('Need an irf or point spread function bintable')
        # get data matrix
        self.get_data_matrices()

    def columns(self):
        """Get FITS extension columns.
        
        Parameter
        ---------
        
        Return
        ------
        columns : list
            FITS extension colums
        """
        return self.psf_data.columns

    def data(self):
        """Get FITS extension data.
        
        Parameter
        ---------

        Return
        ------
        data : ndarray
            FITS extension data
        """
        return self.psf_data.data

    def get_data_matrices(self):
        """Get point spread function data matrix, energy[LO,HI] and theta[LO,HI].

        Parameter
        ---------

        Return
        ------
        psf_matrix : ndarray
            effective area FITS extension matrix
        energies : ndarray
            energy FITS extension subset
        thetas : ndarray
            thetas FITS extension subset
        """
        if self.energies is None or self.thetas is None or self.psf_matrix is None:
            data = self.data()
        if self.energies is None:
            self.energies = np.column_stack((data.field('ENERG_LO')[0], data.field('ENERG_HI')[0]))
        if self.thetas is None:
            self.thetas = np.column_stack((data.field('THETA_LO')[0], data.field('THETA_HI')[0]))
        
        theta_len    = len(self.thetas)
        energies_len = len(self.energies)
        if self.psf_matrix is None:
            fmts = ['f8'] * len(self.fields)
            self.psf_matrix = np.zeros((theta_len, energies_len), dtype={'names': self.fields, 'formats': tuple(fmts) })

            for f in self.fields:
                self.psf_matrix[f] = data.field(f)[0]
        return self.psf_matrix, self.energies, self.thetas

    # maybe a better name: eval_region_psf_rate() or eval_region_rate()
    def eval_region_flux_rate(self, region, pointing, energy):
        """Compute the flux rate in a specific region. This method assumes source region has source as center. The value is the integration from the region center to the edge. If the region contains the full psf (inside 5σ) it returns unit value.

                        1             1    x - μ
        gaussian = ---------- exp( - --- (-------)² )
                   sqrt(2πσ²)         2      σ

        prefactor = 1.0 / (2.0 * np.pi * (sigma_1^2 + ampl_2 * sigma_2^2 + ampl_3 * sigma_3^2))

        Note: 
        - scale value ~= 1.0 / (2.0 * np.pi * (sigma_1 + ampl_2 * sigma_2 + ampl_3 * sigma_3))

        Parameter
        ---------
        region : dict
            region dictionary with ra, dec, rad in degrees
        pointing : dict
            pointing dictionary with ra, dec in degrees

        Return
        ------
        val : float
            integrated flux rate
        """
        region_center = utils.get_skycoord(region)
        region_radius = utils.get_angle(region['rad'])
        pnt_center    = utils.get_skycoord(pointing)
        theta = pnt_center.separation(region_center)

        delta_max = self.get_psf_delta_max(theta, energy)
        if delta_max <= region_radius.degree:
            return (1.0, 0.0)

        sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3 = self.get_psf_values(theta, energy)
        sigmas2_rad = [ np.deg2rad(s)**2 for s in [sigma_1, sigma_2, sigma_3] ]
        prefactor_rad = 1.0 / (2.0 * np.pi * sigmas2_rad[0] + ampl_2 * sigmas2_rad[1] + ampl_3 * sigmas2_rad[2])

        def psf_value(delta_rad):
            """Compute the point spread function value.

            Parameter
            ---------
            delta_rad : float
                delta radius for integration

            Return
            ------
            val : float
                value of the psf
            """
            d2 = delta_rad**2
            numerator  = np.exp( -1/2 * d2 / sigmas2_rad[0] )
            numerator += np.exp( -1/2 * d2 / sigmas2_rad[1] ) * ampl_2 if sigma_2 > 0 else 0
            numerator += np.exp( -1/2 * d2 / sigmas2_rad[2] ) * ampl_3 if sigma_3 > 0 else 0
            return prefactor_rad * numerator

        # integration to 0 to rad on region circumference of psf value
        crf_psf_fn = lambda delta: psf_value(delta) * 2.0 * np.pi * np.sin(delta)
        return integrate.quad(crf_psf_fn, 0, np.deg2rad(region_radius.degree))

    def get_psf_engine(self, region, pointing, energy):
        """Get the psf engine. The engine function can elaborate the psf rate given starting and stop angle [rad]. Each engine depends by theta (between source region and pointing, and energy.

        Note: 
        - delta_max is not the best things to implement in this context because the integration distance is variabile and there is no 'official' source region radius. Furthermore, we could implement the delta_max limit in _integrate_psf considering as radius the difference between integration params and the (eventually) region radius passed as input.

        Parameter
        ---------
        region : dict
            region dictionary with ra, dec, rad in degrees
        pointing : dict
            pointing dictionary with ra, dec in degrees
        energy : float
            energy in [TeV]

        Return
        ------
        val :  float
            integrate point spread function value
        """
        region_center = utils.get_skycoord(region)
        region_radius = utils.get_angle(region['rad'])
        pnt_center    = utils.get_skycoord(pointing)
        theta = pnt_center.separation(region_center)

        sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3 = self.get_psf_values(theta, energy)
        sigmas2_rad = [ np.deg2rad(s)**2 for s in [sigma_1, sigma_2, sigma_3] ]
        prefactor_rad = 1.0 / (2.0 * np.pi * sigmas2_rad[0] + ampl_2 * sigmas2_rad[1] + ampl_3 * sigmas2_rad[2])

        def _integrate_psf(start_rad, stop_rad):
            """Integrate the point spread function.
            
            Parameter
            ---------
            start_rad : float
                center of the region in radians
            stop_rad : fload
                edge of region in radians

            Return
            ------
            val : float
                integrated point spread function value
            """
            if start_rad < 0:
                raise Exception('The starting angle [rad] must be positive')

            def psf_value(delta_rad):
                """Compute the point spread function value.

                Parameter
                ---------
                delta_rad : float
                    delta radius for integration

                Return
                ------
                val : float
                    value of the psf
                """
                d2 = delta_rad**2
                numerator  = np.exp( -1/2 * d2 / sigmas2_rad[0] )
                numerator += np.exp( -1/2 * d2 / sigmas2_rad[1] ) * ampl_2 if sigma_2 > 0 else 0
                numerator += np.exp( -1/2 * d2 / sigmas2_rad[2] ) * ampl_3 if sigma_3 > 0 else 0
                return prefactor_rad * numerator

            # integration to start_rad to stop_rad on region circumference of psf value
            crf_psf_fn = lambda delta: psf_value(delta) * 2.0 * np.pi * np.sin(delta)
            return integrate.quad(crf_psf_fn, start_rad, stop_rad)
        return _integrate_psf

    def get_psf_delta_max(self, offset, energy):
        """Compute maximum delta value [deg] for theta and energy, which is 5*sigma value. 

        Parameter
        ---------
        offset : astropy Angle
            theta angle between center region and pointing
        energy : array
            energy values in [TeV]

        Return
        ------
        val : float
            maximum delta sigma value
        """
        sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3 = self.get_psf_values(offset, energy)
        sigma = sigma_1
        if sigma_2 > sigma:
            sigma = sigma_2
        if sigma_3 > sigma:
            sigma = sigma_3
        return 5.0 * sigma

    def get_psf_values(self, offset, energy):
        """Get the psf data array (sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3). This method returns plain value. No interpolation.

        Parameter
        ---------
        offset : astropy Angle
            offset angle from pointing
        energy : array
            energy values in [TeV]

        Return
        ------
        psf : ndarray
            point spread function ndarray data
        """
        offset_angle = utils.get_angle(offset)
        psf_matrix, energy_bins, theta_bins = self.get_data_matrices()

        theta_index = None
        energy_index = None
        for i, tb in enumerate(theta_bins):
            if offset_angle.degree >= tb[0] and offset_angle.degree < tb[1]:
                theta_index = i
                break
        for j, en in enumerate(energy_bins):
            if energy >= en[0] and energy < en[1]:
                energy_index = j
                break
        if theta_index is None:
            raise Exception('Theta offset is out of range ({})'.format(offset))
        if energy_index is None:
            raise Exception('Energy is out of range ({})'.format(energy))
        return psf_matrix[theta_index, energy_index]

    # this interpolated is a test
    def get_psf_1d_log(self, offset, energy):
        """Compute psf data array (sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3) via 1d interpolation on energy.

        Parameter
        ---------
        offset : astropy Angle
        energy : array 
            energy values in [TeV]

        Return
        ------
        psf : array
            values of interpolated psf
        """
        offset_angle = utils.get_angle(offset)
        psf_matrix, energy_bins, theta_bins = self.get_data_matrices()

        theta_index = None
        for i, tb in enumerate(theta_bins):
            if offset_angle.degree >= tb[0] and offset_angle.degree < tb[1]:
                theta_index = i
                break
        if theta_index is None:
            raise Exception('Theta offset is out of range ({})'.format(offset))

        # energy interpolation
        energy_mid = [ (e[0]+e[1])/2 for e in np.log10(energy_bins) ]
        y_values = []
        for f in self.fields:
            y_values.append(psf_matrix[theta_index][f])
        interp_fn = interpolate.interp1d(x = energy_mid, y = y_values)
        return tuple(interp_fn(np.log10(energy)))

def aeff_eval(args, src, pnt):
    '''Compute the region response of IRF.
    
    Parameter
    ---------
    src : dict
        on region coordinates with ra, dec in degrees
    pnt : dict
        pointing coordinates with ra, dec in degrees

    Return
    ------
    on_reg_aeff : float
        effective area response in [cm2] for on region
    '''
    if not(args.emin and args.emax and args.pixel_size and args.power_law_index):
        raise Exception('need energy min and max, a pixel size to eval the flux')

    aeff = EffectiveArea(irf_filename=args.irf_file)
    # these IRFs return value in m², so we need convert
    # the source data struct need a 'rad'
    on_reg_aeff = aeff.weighted_value_for_region(src, pnt, [args.emin, args.emax], args.pixel_size, args.power_law_index) * 1e4 # cm2
    return on_reg_aeff
