# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Simone Tampieri <simone.tampieri@inaf.it>
# *******************************************************************************

from astropy.io import fits
import numpy as np
from scipy import interpolate, integrate
from rtasci.aph import utils 
from astropy.coordinates import SkyCoord
import math

class IRF:
    def __init__(self, filename):
        self.filename = filename
        self.hdul = fits.open(self.filename)

    def get_extension(self, name):
        return self.hdul[name]

    def get_eff_area(self):
        return self.get_extension('EFFECTIVE AREA')

    def get_psf_data(self):
        return self.get_extension('POINT SPREAD FUNCTION')

class EffectiveArea:
    def __init__(self, irf_filename=None, eff_area_bintable=None):
        self.irf_filename = None
        self.eff_area = None

        # a sort of cache...
        self.energies = None
        self.thetas = None
        self.aeff_matrix = None

        if irf_filename is not None:
            self.irf_filename = irf_filename
            irf = IRF(self.irf_filename)
            self.eff_area = irf.get_eff_area()
        elif eff_area_bintable is not None:
            self.eff_area = eff_area_bintable

        if self.eff_area is None:
            raise Exception('Need an irf or effective area bintable')

        self.get_data_matrices()

    def columns(self):
        return self.eff_area.columns

    def data(self):
        return self.eff_area.data

    def get_data_matrices(self):
        """ returns aeff data matrix, energy[LO,HI] and theta[LO,HI]
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
        """
        return effective area in [m²]

        Parameters
          offset: Angle
          energy: [TeV]

        This method does 1D interpolation on energy range, managed as log10.
        Theta offset is not interpolated.
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
        """
        return effective area in [m²]

        Parameters
          offset: Angle
          energy: [TeV]

        This method does 2D interpolation and manage energy as log10.
        """
        offset_angle = utils.get_angle(input_offset)
        aeff_matrix, energy_bins, theta_bins = self.get_data_matrices()

        # energy interpolation
        theta_mid  = [ (t[0]+t[1])/2 for t in theta_bins ]
        energy_mid = [ (e[0]+e[1])/2 for e in np.log10(energy_bins) ]
        energy_fn  = interpolate.interp2d(x = energy_mid, y = theta_mid, z = aeff_matrix)
        return energy_fn(np.log10(input_energy), input_offset)[0]

    def weighted_value_for_region(self, *args):
        return self.weighted_aeff_flat_psf_w_powerlaw(*args)

    # FIXME the idea is integrate psf gradually but the offset from center and from pointing works differently
    # we need to associate the PSF degradation from source region center
    # and the Aeff degradation form the pointing
    # def weighted_aeff_psf_w_powerlaw(self, region, pointing, input_energies, pixel_size=0.05, e_index=-2.4):
    #    ...

    # this method use the psf value plain as it easy for energy weight. doesn't
    # consider the spatial component of the PSF
    def weighted_aeff_flat_psf_w_powerlaw(self, region, pointing, input_energies, pixel_size=0.05, e_index=-2.4):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energies: a couple of values in TeV (ex: [ 0.025, 1.0 ])
          pixel_size: a value in degree (default: 0.05)
          e_index: is the powerlaw index (default: -2.4)
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

    # this method use an energy range to evaluate the aeff.
    # The energy range is binned and weighted with a powerlaw with index = e_index.
    # Each pixel has a radial weigth and a specific column of energies weight.
    # Lower energies have more weigth than higher.
    # if energy range is small, the effect is trascurable - similar to weighted_value_for_region_single_energy method.
    def weighted_value_for_region_w_powerlaw(self, region, pointing, input_energies, pixel_size=0.05, e_index=-2.4):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energies: a couple of values in TeV (ex: [ 0.025, 1.0 ])
          pixel_size: a value in degree (default: 0.05)
          e_index: is the powerlaw index (default: -2.4)
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

    # this method use an energy range to evaluate the aeff. The energy range is
    # binned and every matrix cube (pixel distance * energy bin) have the same
    # weight.
    # no powerlaw here
    def weighted_value_for_region_no_powerlaw(self, region, pointing, input_energies, pixel_size=0.05):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energies: a couple of values in TeV (ex: [ 0.025, 1.0 ])
          pixel_size: a value in degree (default: 0.05)
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

    # Deprecated 2019-12-11
    # this method is old. just a plain output from an array of points and ONE
    # energy (usually the middle point of the energy range)
    def weighted_value_for_region_single_energy(self, region, pointing, energy, pixel_size=0.05):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energy: a value in TeV (ex: 0.025, 1.0, 150.0)
          pixel_size: a value in degree (default: 0.05)
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
    def __init__(self, irf_filename=None, psf_bintable=None):
        self.irf_filename = None
        self.psf_data = None
        self.fields = ('SIGMA_1', 'SIGMA_2', 'SIGMA_3', 'SCALE', 'AMPL_2', 'AMPL_3')

        # a sort of cache...
        self.energies = None
        self.thetas = None
        self.psf_matrix = None

        if irf_filename is not None:
            self.irf_filename = irf_filename
            irf = IRF(self.irf_filename)
            self.psf_data = irf.get_psf_data()
        elif psf_bintable is not None:
            self.psf_data = psf_bintable

        if self.psf_data is None:
            raise Exception('Need an irf or point spread function bintable')

        self.get_data_matrices()

    def columns(self):
        return self.psf_data.columns

    def data(self):
        return self.psf_data.data

    def get_data_matrices(self):
        """ returns psf data matrix, energy[LO,HI] and theta[LO,HI]
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
        """
        return flux rate in specific source region and integral computation error

        Parameters
            region: source region (ra, dec, rad)
            pointing: pointing direction (ra, dec)

        This method assumes source region has source as center. Is not general
        enough, but at the moment is ok.
        The value is the integration from the region center to the edge. If region
        contains the full psf (inside 5σ) we return 1.0 .

                        1             1    x - μ
        gaussian = ---------- exp( - --- (-------)² )
                   sqrt(2πσ²)         2      σ

        prefactor = 1.0 / (2.0 * np.pi * (sigma_1^2 + ampl_2 * sigma_2^2 + ampl_3 * sigma_3^2))
        Note: 
          scale value ~= 1.0 / (2.0 * np.pi * (sigma_1 + ampl_2 * sigma_2 + ampl_3 * sigma_3))
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
            d2 = delta_rad**2
            numerator  = np.exp( -1/2 * d2 / sigmas2_rad[0] )
            numerator += np.exp( -1/2 * d2 / sigmas2_rad[1] ) * ampl_2 if sigma_2 > 0 else 0
            numerator += np.exp( -1/2 * d2 / sigmas2_rad[2] ) * ampl_3 if sigma_3 > 0 else 0
            return prefactor_rad * numerator

        # integration to 0 to rad on region circumference of psf value
        crf_psf_fn = lambda delta: psf_value(delta) * 2.0 * np.pi * np.sin(delta)
        return integrate.quad(crf_psf_fn, 0, np.deg2rad(region_radius.degree))

    def get_psf_engine(self, region, pointing, energy):
        """
        return psf engine. The engine function can elaborate the psf rate given
        starting and stop angle [rad]. Each engine depends by theta (between
        source region and pointing, and energy.

        Parameters
            region: source region (ra, dec, rad)
            pointing: pointing direction (ra, dec)
            energy: energy in TeV
        """
        region_center = utils.get_skycoord(region)
        region_radius = utils.get_angle(region['rad'])
        pnt_center    = utils.get_skycoord(pointing)
        theta = pnt_center.separation(region_center)

        # Note: delta_max is not the best things to implement in this context
        #       'cause the integration distance is variabile and there is no
        #       "official" source region radius.
        #       Furthermore, we could implement the delta_max limit in _integrate_psf
        #       considering as radius the difference between integration params
        #       and the (eventually) region.radius passed as input.
        # delta_max = self.get_psf_delta_max(theta, energy)
        # if delta_max <= region_radius.degree:
        #     return (1.0, 0.0)

        sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3 = self.get_psf_values(theta, energy)
        sigmas2_rad = [ np.deg2rad(s)**2 for s in [sigma_1, sigma_2, sigma_3] ]
        prefactor_rad = 1.0 / (2.0 * np.pi * sigmas2_rad[0] + ampl_2 * sigmas2_rad[1] + ampl_3 * sigmas2_rad[2])

        def _integrate_psf(start_rad, stop_rad):
            # typical params:
            #   start_rad = 0
            #   stop_rad = np.deg2rad(region_radius.degree))
            if start_rad < 0:
                raise Exception('The starting angle [rad] must be positive')
            def psf_value(delta_rad):
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
        """
        return max delta value [deg] for theta and energy.

        Parameters
          offset: Angle
          energy: [TeV]

        max delta value is = 5*sigma value. 
        """
        sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3 = self.get_psf_values(offset, energy)
        sigma = sigma_1
        if sigma_2 > sigma:
            sigma = sigma_2
        if sigma_3 > sigma:
            sigma = sigma_3
        return 5.0 * sigma

    def get_psf_values(self, offset, energy):
        """
        return psf data array (sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3)

        Parameters
          offset: Angle
          energy: [TeV]

        This method returns plain value. No interpolation.
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
        """
        return psf data array (sigma_1, sigma_2, sigma_3, scale, ampl_2, ampl_3)

        Parameters
          offset: Angle
          energy: [TeV]

        This method does 1D interpolation on energy range, managed as log10.
        Theta offset is not interpolated.
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

